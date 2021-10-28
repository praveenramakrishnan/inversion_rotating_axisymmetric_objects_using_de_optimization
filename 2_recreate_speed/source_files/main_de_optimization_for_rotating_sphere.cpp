/*
Author: Praveen Kalarickel Ramakrishnan
Date: 2021-10-27

This program uses De optimization algorithm to get the solution for the 
relative permittivity and the rotating speed of homogelneous sphere.

The interface to th DE optimization algorithm consits of four functions
1. NumberOfParameters()
        Pass in the paramter inputvariable::number_of_unknowns
2. GetConstraints()
        Pass in the vector of unknows
          global_variable::min_constraint[N]
          global_variable::max_constraint[N]
3. termination_condition()
        Depends on the variables:
            inputvariable::max_num_iterations, 
            inputvariable::max_num_iterations_without_improvement,
            inputvariable::max_error_tolerrance_without_improvement.
4. EvaluteCost()
            Needs the output files containing the actual solution.
                Define them as global variables
            Needs a function to calculate the forward problem for test problem.
                Write a function called solveForwardProble(test_inputs)
            Needs the inputs required to solve the forward problem.
                Write a function called errorInSolution(testSolution)

The steps for obtaining the solution are as follows:
1) Read in the input variables

    1.Minumum and maximum values of $\beta$ and $\varepsilonr$ where the solution has to be searched.
*/

#include <functional>
#include "de/DifferentialEvolution.h"
#include <ctime>
#include <string>
#include <fstream>

//Global Variables to be read in as inputs to the program
   namespace input_variable{
      int problem_type;
      double beta_min, beta_max, epsilonr_min, epsilonr_max;
      double max_error_tolerrance_without_improvement;
      int max_num_iterations, max_num_iterations_without_improvement, population_size, num_of_points;
      std::string inputfilename_efield, inputfilename_hfield;
      std::vector<double> real_x_component_efield, imag_x_component_efield, 
                          real_y_component_efield, imag_y_component_efield, 
                          real_z_component_efield, imag_z_component_efield;
      std::vector<double> real_x_component_hfield, imag_x_component_hfield, 
                          real_y_component_hfield, imag_y_component_hfield, 
                          real_z_component_hfield, imag_z_component_hfield;
      double mu_r, epsilonr_true_value, R_sphere, frequency, beta_true_value, theta_i, SNR;
      std::string outputfilename_points;
   }

   namespace global_variable{
      int current_iteration_number = 0, number_of_cost_function_evaluations = 0,
          number_of_unknowns;
      double mean_square_efield, mean_square_hfield;
      std::vector<double> cost_vector, min_constraint, max_constraint;
   }

    //Function to read input data files containing electric and magentic fields
    void read_input_data(std::string inputfilename, int N_points, 
                         std::vector<double>& real_x_component, std::vector<double>& imag_x_component, 
		                 std::vector<double>& real_y_component, std::vector<double>& imag_y_component, 
			             std::vector<double>& real_z_component, std::vector<double>& imag_z_component);

    //Function to calculate the difference between two vectors of same size
    std::vector<double> difference_vector(std::vector<double>& vector1, std::vector<double>& vector2);

    //Function to calculate the sum of squares of components of a vector
    double sum_of_squares(std::vector<double>& vector1);
   
    //Function to calculate sum of squares of 3d complex vector
    double sum_of_squares_3d_complex(std::vector<double>& vector_x_real, std::vector<double>& vector_x_imag, 
		                            std::vector<double>& vector_y_real, std::vector<double>& vector_y_imag,
		                            std::vector<double>& vector_z_real, std::vector<double>& vector_z_imag);
   
   void calculate_forward_solution(double beta_test,  double epsilonr_test,
                                   std::string outputfilename_fwd_problem_efield, 
                                   std::string outputfilename_fwd_problem_hfield);
   double calculate_mean_square_error(std::string outputfilename_fwd_problem_efield, 
                                      std::string outputfilename_fwd_problem_hfield);

   class CostFunction : public de::IOptimizable
   {
      public:
      double EvaluteCost(std::vector<double> inputs) const override
      {
          assert(inputs.size() == NumberOfParameters());

	      global_variable::number_of_cost_function_evaluations += 1;

          //Select according to problem type
          double beta_test, epsilonr_test;
          if(input_variable::problem_type == 0) //beta is the only unknown
          {
            beta_test = inputs[0];
            epsilonr_test = input_variable::epsilonr_true_value;
          }
          else if(input_variable::problem_type == 1) // epsilonr is the only unknown
          {
            beta_test = input_variable::beta_true_value;
            epsilonr_test = inputs[global_variable::number_of_unknowns -1];
          }
          else if(input_variable::problem_type == 2) // Both beta and epsilonr are unknowns
          {
            beta_test = inputs[0];
            epsilonr_test = inputs[global_variable::number_of_unknowns -1];
          }
	      //Calculate the fields for given beta value
          std::string outputfilename_fwd_problem_efield = "test_solution_efield";
	      std::string outputfilename_fwd_problem_hfield = "test_solution_hfield";
          calculate_forward_solution(beta_test, epsilonr_test, 
                                     outputfilename_fwd_problem_efield, outputfilename_fwd_problem_hfield);
	      
	      //Calculate the mean square difference between forward solution and the measured values
          return calculate_mean_square_error(outputfilename_fwd_problem_efield, outputfilename_fwd_problem_hfield);
	      //system("rm test_solution_efield test_solution_hfield");
      }
	

    unsigned int NumberOfParameters() const override
    {
        return global_variable::number_of_unknowns;
    }

    std::vector<Constraints> GetConstraints() const override
    {
        std::vector<Constraints> constr(NumberOfParameters());
        
        for(int iloop = 0; iloop < NumberOfParameters(); iloop++)
        {
            constr[iloop] = Constraints(global_variable::min_constraint[iloop], global_variable::max_constraint[iloop], true);
        }
       return constr;
    }

   };

    bool termination_condition(const de::DifferentialEvolution& de)
    {
    	global_variable::cost_vector.push_back(de.m_get_minCost());
    	global_variable::current_iteration_number += 1;
    
    	int iter_num = global_variable::current_iteration_number-1;
    	int prev_iter_num = global_variable::current_iteration_number - input_variable::max_num_iterations_without_improvement;
    
    	if(prev_iter_num >= 0)
    	{
             //if(iter_num > input_variable:;max_num_iterations)
             //{
             //    return true;
             //}
    
             double error_improvement =  (global_variable::cost_vector[prev_iter_num] 
                                             - global_variable::cost_vector[iter_num])
                                             / global_variable::cost_vector[iter_num];
    
    	     if(error_improvement < input_variable::max_error_tolerrance_without_improvement) 
    	     {
    		        return true;
    	     }
    	}
    
    	return false;
    }

    void print_output_file();

int main()
{
     //Read in the input variables
     std::cin >> input_variable::problem_type;
     if(input_variable::problem_type == 0)
     {
        std::cin  >> input_variable::beta_min >> input_variable::beta_max;  
        global_variable::number_of_unknowns = 1;
     }
     else if(input_variable::problem_type == 1)
     {
        std::cin  >> input_variable::epsilonr_min >> input_variable::epsilonr_max;
        global_variable::number_of_unknowns = 1;
     }
     else if(input_variable::problem_type == 2)
     {
        std::cin  >> input_variable::beta_min >> input_variable::beta_max;  
        std::cin  >> input_variable::epsilonr_min >> input_variable::epsilonr_max;
        global_variable::number_of_unknowns = 2;
     }
	 std::cin >> input_variable::max_num_iterations >> input_variable::max_num_iterations_without_improvement 
	          >> input_variable::max_error_tolerrance_without_improvement >> input_variable::population_size  
	          >> input_variable::num_of_points >> input_variable::inputfilename_efield >> input_variable::inputfilename_hfield;

     std::cin >> input_variable::mu_r >> input_variable::epsilonr_true_value >> input_variable::R_sphere 
              >> input_variable::frequency >> input_variable::beta_true_value 
              >> input_variable::theta_i >> input_variable::outputfilename_points >> input_variable::SNR;

    //Read inputfile containing electric field
    read_input_data(input_variable::inputfilename_efield, input_variable::num_of_points, 
                    input_variable::real_x_component_efield, input_variable::imag_x_component_efield, 
		            input_variable::real_y_component_efield, input_variable::imag_y_component_efield, 
                    input_variable::real_z_component_efield, input_variable::imag_z_component_efield);


    //Read inputfile containing magnetic field
    read_input_data(input_variable::inputfilename_hfield, input_variable::num_of_points, 
                    input_variable::real_x_component_hfield, input_variable::imag_x_component_hfield, 
		            input_variable::real_y_component_hfield, input_variable::imag_y_component_hfield, 
                    input_variable::real_z_component_hfield, input_variable::imag_z_component_hfield);


    //Calculate measn square value of electric field
    global_variable::mean_square_efield = sum_of_squares_3d_complex(
                                              input_variable::real_x_component_efield, input_variable::imag_x_component_efield, 
		                                      input_variable::real_y_component_efield, input_variable::imag_y_component_efield,
		                                      input_variable::real_z_component_efield, input_variable::imag_z_component_efield);
    assert(global_variable::mean_square_efield != 0);

    //Calculate mean square value of magnetic field
    global_variable::mean_square_hfield = sum_of_squares_3d_complex(
                                               input_variable::real_x_component_hfield, input_variable::imag_x_component_hfield, 
		                                       input_variable::real_y_component_hfield, input_variable::imag_y_component_hfield,
		                                       input_variable::real_z_component_hfield, input_variable::imag_z_component_hfield);
    assert(global_variable::mean_square_hfield != 0);

    //make the constraint vector
    global_variable::min_constraint.resize(global_variable::number_of_unknowns);
    global_variable::max_constraint.resize(global_variable::number_of_unknowns);
    if(input_variable::problem_type == 0 || input_variable::problem_type == 2)
     {
        global_variable::min_constraint[0] = input_variable::beta_min;
        global_variable::max_constraint[0] = input_variable::beta_max;
     }
     if(input_variable::problem_type == 1 || input_variable::problem_type == 2)
     {
        global_variable::min_constraint[global_variable::number_of_unknowns -1] = input_variable::epsilonr_min;
        global_variable::max_constraint[global_variable::number_of_unknowns -1] = input_variable::epsilonr_max;
     }
    
    //Print the input parameters to the output file

    print_output_file();

    //Running differential evolution algorithm
    
    CostFunction cost;

    de::DifferentialEvolution de(cost, input_variable::population_size, std::time(nullptr), true, nullptr, termination_condition);
    
    de.Optimize(input_variable::max_num_iterations, true);

    //print output solution to outputfile
    
    std::cout << "Number of cost function evaluations = " << global_variable::number_of_cost_function_evaluations << std::endl;

    if(input_variable::problem_type == 0 || input_variable::problem_type == 2)
    {
        double error_in_beta = std::abs(de.m_get_bestAgent()[0] - input_variable::beta_true_value)/input_variable::beta_true_value;

        std::cout << "solution for beta = " << de.m_get_bestAgent()[0] << std::endl;

        std::cout << "error in beta = " << error_in_beta << std::endl;
    }

    if(input_variable::problem_type == 1 || input_variable::problem_type == 2)
    {
        double error_in_epsilonr = std::abs(de.m_get_bestAgent()[global_variable::number_of_unknowns-1] 
                                 - input_variable::epsilonr_true_value)
                                 /input_variable::epsilonr_true_value;

        std::cout << "solution for epsilonr  = " << de.m_get_bestAgent()[global_variable::number_of_unknowns-1] << std::endl;

        std::cout << "error in epsilonr = " << error_in_epsilonr << std::endl;
    }

    return 0;
}


   //Function to read input data files containing electric and magentic fields
   void read_input_data(std::string inputfilename, int N_points, 
                         std::vector<double>& real_x_component, std::vector<double>& imag_x_component, 
		                 std::vector<double>& real_y_component, std::vector<double>& imag_y_component, 
			             std::vector<double>& real_z_component, std::vector<double>& imag_z_component)
    {
        std::ifstream inputfile;
        inputfile.open(inputfilename.c_str());
        assert(inputfile.is_open());

        double dummy_double; 
	    std::vector<double> x_coord(N_points), y_coord(N_points), z_coord(N_points); 
        real_x_component.resize(N_points), imag_x_component.resize(N_points),
        real_y_component.resize(N_points), imag_y_component.resize(N_points), 
        real_z_component.resize(N_points), imag_z_component.resize(N_points);


        for(int i_loop=0; i_loop < N_points; i_loop++)
        {
           inputfile >> x_coord[i_loop] >> y_coord[i_loop] >> z_coord[i_loop] 
           	         >> real_x_component[i_loop] >> imag_x_component[i_loop] >> dummy_double >> dummy_double
                     >> real_y_component[i_loop] >> imag_y_component[i_loop] >> dummy_double >> dummy_double
             	     >> real_z_component[i_loop] >> imag_z_component[i_loop] >> dummy_double >> dummy_double;
        }
        inputfile.close();

    }

    //Function to calculate the difference between two vectors of same size
    std::vector<double> difference_vector(std::vector<double>& vector1, std::vector<double>& vector2)
    {
        int N_points = vector1.size();
        assert(N_points == vector2.size());
        std::vector<double> vector_difference(N_points);

        for(int i_loop=0; i_loop < N_points; i_loop++){
              vector_difference[i_loop] = vector1[i_loop] - vector2[i_loop];
        }

        return vector_difference;
    }

    //Function to calculate the sum of squares of components of a vector
    double sum_of_squares(std::vector<double>& vector1)
    {
    	int N_points = vector1.size();
    	double square_sum=0;
    
        for(int i_loop=0; i_loop < N_points; i_loop++){
    		   square_sum += vector1[i_loop]*vector1[i_loop];
    	}
    
    	return square_sum;
    }

    //Function to calculate sum of squares of 3d complex vector
    double sum_of_squares_3d_complex(std::vector<double>& vector_x_real, std::vector<double>& vector_x_imag, 
		                            std::vector<double>& vector_y_real, std::vector<double>& vector_y_imag,
		                            std::vector<double>& vector_z_real, std::vector<double>& vector_z_imag)
    {
	       return sum_of_squares(vector_x_real) + sum_of_squares(vector_x_imag) 
	             + sum_of_squares(vector_y_real) + sum_of_squares(vector_y_imag) 
	             + sum_of_squares(vector_z_real) + sum_of_squares(vector_z_imag);
    }

   void calculate_forward_solution(double beta_test, double epsilonr_test,
                                   std::string outputfilename_fwd_problem_efield, 
                                   std::string outputfilename_fwd_problem_hfield)
   {
          std::string infilename_fwd_problem = "infile_dezutter";
	      
	      std::ofstream inputfile_fwd_problem;
          inputfile_fwd_problem.open(infilename_fwd_problem.c_str());
          assert(inputfile_fwd_problem.is_open());

          inputfile_fwd_problem <<  std::scientific << std::setprecision(16);
	      inputfile_fwd_problem <<  input_variable::mu_r << std::endl;
	      //inputfile_fwd_problem <<  input_variable::epsilonr_true_value << std::endl;
	      inputfile_fwd_problem <<  epsilonr_test << std::endl;
	      inputfile_fwd_problem <<  input_variable::R_sphere << std::endl;
	      inputfile_fwd_problem <<  input_variable::frequency << std::endl;
	      inputfile_fwd_problem <<  beta_test << std::endl;
	      inputfile_fwd_problem <<  input_variable::theta_i<< std::endl;
	      inputfile_fwd_problem <<  input_variable::outputfilename_points << std::endl;
	      inputfile_fwd_problem <<  outputfilename_fwd_problem_efield << std::endl;
	      inputfile_fwd_problem <<  outputfilename_fwd_problem_hfield << std::endl;

          inputfile_fwd_problem.close();

	      system("time ./main_rotating_sphere_dezutter_analytic_solution.o < infile_dezutter; rm infile_dezutter");
   }

   double calculate_mean_square_error(std::string outputfilename_fwd_problem_efield, 
                                      std::string outputfilename_fwd_problem_hfield)
   {
          //Read in the calculated test solutions
	      std::vector<double> real_x_component_efield_test, imag_x_component_efield_test, 
		                      real_y_component_efield_test, imag_y_component_efield_test,
                              real_z_component_efield_test, imag_z_component_efield_test;

          read_input_data(outputfilename_fwd_problem_efield, input_variable::num_of_points, 
                          real_x_component_efield_test, imag_x_component_efield_test, 
		                  real_y_component_efield_test, imag_y_component_efield_test, 
                          real_z_component_efield_test, imag_z_component_efield_test);

	      std::vector<double> real_x_component_hfield_test, imag_x_component_hfield_test, 
		                      real_y_component_hfield_test, imag_y_component_hfield_test,
                              real_z_component_hfield_test, imag_z_component_hfield_test;

          read_input_data(outputfilename_fwd_problem_hfield, input_variable::num_of_points, 
                          real_x_component_hfield_test, imag_x_component_hfield_test, 
		                  real_y_component_hfield_test, imag_y_component_hfield_test, 
                          real_z_component_hfield_test, imag_z_component_hfield_test);
	
	     //Find the mean square difference between calculated fields and "measured" fields
	     std::vector<double> difference_vector_real_x_component_efield = 
                                     difference_vector(real_x_component_efield_test, input_variable::real_x_component_efield);
	     std::vector<double> difference_vector_imag_x_component_efield = 
                                     difference_vector(imag_x_component_efield_test, input_variable::imag_x_component_efield);
	     std::vector<double> difference_vector_real_y_component_efield = 
                                     difference_vector(real_y_component_efield_test, input_variable::real_y_component_efield);
	     std::vector<double> difference_vector_imag_y_component_efield = 
                                     difference_vector(imag_y_component_efield_test, input_variable::imag_y_component_efield);
	     std::vector<double> difference_vector_real_z_component_efield = 
                                     difference_vector(real_z_component_efield_test, input_variable::real_z_component_efield);
	     std::vector<double> difference_vector_imag_z_component_efield = 
                                     difference_vector(imag_z_component_efield_test, input_variable::imag_z_component_efield);

	     double mean_square_efield_error = (sum_of_squares_3d_complex(
                                               difference_vector_real_x_component_efield, difference_vector_imag_x_component_efield, 
                                               difference_vector_real_y_component_efield, difference_vector_imag_y_component_efield, 
                                               difference_vector_real_z_component_efield, difference_vector_imag_z_component_efield) ) 
	     	                                /global_variable::mean_square_efield;

	     std::vector<double> difference_vector_real_x_component_hfield = 
                                     difference_vector(real_x_component_hfield_test, input_variable::real_x_component_hfield);
	     std::vector<double> difference_vector_imag_x_component_hfield = 
                                     difference_vector(imag_x_component_hfield_test, input_variable::imag_x_component_hfield);
	     std::vector<double> difference_vector_real_y_component_hfield = 
                                     difference_vector(real_y_component_hfield_test, input_variable::real_y_component_hfield);
	     std::vector<double> difference_vector_imag_y_component_hfield = 
                                     difference_vector(imag_y_component_hfield_test, input_variable::imag_y_component_hfield);
	     std::vector<double> difference_vector_real_z_component_hfield = 
                                     difference_vector(real_z_component_hfield_test, input_variable::real_z_component_hfield);
	     std::vector<double> difference_vector_imag_z_component_hfield = 
                                     difference_vector(imag_z_component_hfield_test, input_variable::imag_z_component_hfield);

	     double mean_square_hfield_error = (sum_of_squares_3d_complex(
                                               difference_vector_real_x_component_hfield, difference_vector_imag_x_component_hfield, 
                                               difference_vector_real_y_component_hfield, difference_vector_imag_y_component_hfield, 
                                               difference_vector_real_z_component_hfield, difference_vector_imag_z_component_hfield) ) 
	     	                                /global_variable::mean_square_hfield;


        //std::string system_command = "rm test_solution_efield test_solution_hfield";
	    //system(system_command.c_str());
	    system("rm test_solution_efield test_solution_hfield");
	
        return mean_square_efield_error + mean_square_hfield_error; 
   }

    void print_output_file()
    {
        std::cout << "number of evaluation points = " << input_variable::num_of_points << std::endl;
        std::cout << "permitivity = " << input_variable::epsilonr_true_value << std::endl;
        std::cout << "permeability = " << input_variable::mu_r << std::endl;
        std::cout << "Radius of sphere = " << input_variable::R_sphere << std::endl;
        std::cout << "True beta value = " << input_variable::beta_true_value << std::endl;
        std::cout << "Angle of incidence = " << input_variable::theta_i << std::endl;
        if(input_variable::problem_type == 0 || input_variable::problem_type == 2)
        {
            std::cout << "beta minimum = " << input_variable::beta_min << std::endl;
            std::cout << "beta maximum = " << input_variable::beta_max << std::endl;
        }
        if(input_variable::problem_type == 0 || input_variable::problem_type == 2)
        {
            std::cout << "epsilonr minimum = " << input_variable::epsilonr_min << std::endl;
            std::cout << "epsilonr maximum = " << input_variable::epsilonr_max << std::endl;
        }
        std::cout << "maximum number of iterations = " << input_variable::max_num_iterations << std::endl;
        std::cout << "maximum number of iterations without improvement = " 
                  << input_variable::max_num_iterations_without_improvement << std::endl;
        std::cout << "error tolerance without improvement = " << input_variable::max_error_tolerrance_without_improvement << std::endl;
        std::cout << "population size = " << input_variable::population_size << std::endl;
        std::cout << "SNR = " << input_variable::SNR << std::endl;
    }
