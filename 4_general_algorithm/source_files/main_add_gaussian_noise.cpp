/*
DATE: 05/02/2021
AUTHOR: Praveen Kalarickel Ramakrishnan

This program adds guassian noise of specified SNR to the data of electric or magnetic fields.
*/
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>
#include <random>
#include <iomanip>

using namespace std;
typedef complex<double> dcomplex;
#define j dcomplex(0.0,1.0)
const double pi = 3.1415926535897932;

int number_of_lines_of_file(string inputfilename);
double mean_value_3d_vector(int N_points, vector<double> magnitude_x, vector<double> magnitude_y, vector<double> magnitude_z);

int main()
{
   string inputfilename_field, outputfilename_field;
   double SNR_dB_value, signal_field_level, dummy_double;

   //Read the input data
   cin >> inputfilename_field;

   ifstream inputfile;
   inputfile.open(inputfilename_field.c_str());
   assert(inputfile.is_open());

   int N_points;
   N_points = number_of_lines_of_file(inputfilename_field);

   vector<double> x_coord(N_points), y_coord(N_points), z_coord(N_points), 
                  real_x_component_field(N_points), imag_x_component_field(N_points), mag_x_component_field(N_points),
	              real_y_component_field(N_points), imag_y_component_field(N_points), mag_y_component_field(N_points),
                  real_z_component_field(N_points), imag_z_component_field(N_points), mag_z_component_field(N_points);


   for(int i_loop=0; i_loop < N_points; i_loop++)
   {
      inputfile >> x_coord[i_loop] >> y_coord[i_loop] >> z_coord[i_loop] 
                >> real_x_component_field[i_loop] >> imag_x_component_field[i_loop] >> mag_x_component_field[i_loop] >> dummy_double
	            >> real_y_component_field[i_loop] >> imag_y_component_field[i_loop] >> mag_y_component_field[i_loop] >> dummy_double
		        >> real_z_component_field[i_loop] >> imag_z_component_field[i_loop] >> mag_z_component_field[i_loop] >> dummy_double;
   }
   inputfile.close();

   //Add gaussian noie to input data
   double mean_noise = 0;
   cin >> SNR_dB_value;
   //cin >> signal_field_level;
   signal_field_level = mean_value_3d_vector(N_points, mag_x_component_field, mag_y_component_field, mag_z_component_field);
   double SNR = pow(10, SNR_dB_value/20.0);
   double noise_level = signal_field_level/SNR;

   double mean=0, std_dev=noise_level;
   default_random_engine myseed(time(0));
   normal_distribution <double> my_distribution(mean, std_dev);
  
   vector<double>  real_x_component_field_noisy(N_points), imag_x_component_field_noisy(N_points), real_y_component_field_noisy(N_points), imag_y_component_field_noisy(N_points),
	           real_z_component_field_noisy(N_points), imag_z_component_field_noisy(N_points);

   for(int i_loop=0; i_loop < N_points; i_loop++)
   {
      real_x_component_field_noisy[i_loop] = real_x_component_field[i_loop] + my_distribution(myseed)/sqrt(2.0);
      imag_x_component_field_noisy[i_loop] = imag_x_component_field[i_loop] + my_distribution(myseed)/sqrt(2.0);
      real_y_component_field_noisy[i_loop] = real_y_component_field[i_loop] + my_distribution(myseed)/sqrt(2.0);
      imag_y_component_field_noisy[i_loop] = imag_y_component_field[i_loop] + my_distribution(myseed)/sqrt(2.0);
      real_z_component_field_noisy[i_loop] = real_z_component_field[i_loop] + my_distribution(myseed)/sqrt(2.0);
      imag_z_component_field_noisy[i_loop] = imag_z_component_field[i_loop] + my_distribution(myseed)/sqrt(2.0);
   }

   //Print the output file
   cin >> outputfilename_field;
   ofstream outputfile;
   outputfile.open(outputfilename_field.c_str());
   assert(outputfile.is_open());

   outputfile << scientific << setprecision(16);
   for(int i_loop=0; i_loop < N_points; i_loop++){
      outputfile << x_coord[i_loop] << "\t" << y_coord[i_loop] << "\t" << z_coord[i_loop] <<  "\t" 
	         << real_x_component_field_noisy[i_loop] << "\t" << imag_x_component_field_noisy[i_loop] << "\t" 
		 << sqrt(real_x_component_field_noisy[i_loop]*real_x_component_field_noisy[i_loop] + imag_x_component_field_noisy[i_loop]*imag_x_component_field_noisy[i_loop]) << "\t"
		 << atan2(imag_x_component_field_noisy[i_loop], real_x_component_field_noisy[i_loop]) << "\t"
	         << real_y_component_field_noisy[i_loop] << "\t" << imag_y_component_field_noisy[i_loop] << "\t" 
		 << sqrt(real_y_component_field_noisy[i_loop]*real_y_component_field_noisy[i_loop] + imag_y_component_field_noisy[i_loop]*imag_y_component_field_noisy[i_loop]) << "\t"
		 << atan2(imag_y_component_field_noisy[i_loop], real_y_component_field_noisy[i_loop]) << "\t"
	         << real_z_component_field_noisy[i_loop] << "\t" << imag_z_component_field_noisy[i_loop] << "\t" 
		 << sqrt(real_z_component_field_noisy[i_loop]*real_z_component_field_noisy[i_loop] + imag_z_component_field_noisy[i_loop]*imag_z_component_field_noisy[i_loop]) << "\t"
		 << atan2(imag_z_component_field_noisy[i_loop], real_z_component_field_noisy[i_loop]) << endl;

   }
   outputfile.close();


   return 0;
}

int number_of_lines_of_file(string inputfilename)
{
   int N_lines = 0;
   string dummy_string;

   ifstream inputfile;
   inputfile.open(inputfilename.c_str());
   assert(inputfile.is_open());
   while(getline(inputfile, dummy_string))
	   ++N_lines;
   inputfile.close();
   return N_lines;

}

double mean_value_3d_vector(int N_points, vector<double> magnitude_x, vector<double> magnitude_y, vector<double> magnitude_z)
{
   double mean = 0;

   for(int iloop=0; iloop < N_points; iloop++)
   {

       mean += magnitude_x[iloop]*magnitude_x[iloop] 
             + magnitude_y[iloop]*magnitude_y[iloop] 
             + magnitude_z[iloop]*magnitude_z[iloop];
   }
   mean = sqrt(mean)/N_points;

   return mean;
}
