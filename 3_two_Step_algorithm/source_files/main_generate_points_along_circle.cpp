/* Date: 05/02/2021
   Author: Praveen K R
This program generates points on a circle

Inputs:
Number of points
Plane along which points are printed: 1 for xy, 2 for yz, 3 for zx
Radius
Name of the output file

Output file format:
First line: number of points
Subsequent lines: x, y, z
  */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>

using namespace std;
typedef complex<double> dcomplex;
#define j dcomplex(0.0,1.0)
const double pi = 3.1415926535897932;
const double C0 = 299792458.0;//speed of light in vacuum
const double MU0 = 4.0*pi*1.0e-7; 
const double EPS0 = 1.0/(MU0*C0*C0);

void calculate_points(int N_points, double radius, vector<double>& coord1, vector<double>& coord2, vector<double>& coord3);

int main()
{

   //Read input 
   int N_points, plane_of_circle;
   double radius;

   cin >> N_points;
   cin >> plane_of_circle;
   cin >> radius;

   //Calculate the points 
   vector<double> x_coord(N_points), y_coord(N_points), z_coord(N_points);    

   if( plane_of_circle == 1)
   {
      calculate_points(N_points, radius, x_coord, y_coord, z_coord);
   }
   else if(plane_of_circle == 2)
   {
      calculate_points(N_points, radius, y_coord, z_coord, x_coord);
   }
   else if(plane_of_circle == 3)
   {
      calculate_points(N_points, radius, z_coord, x_coord, y_coord);
   }
   else
   {
      assert(plane_of_circle == 1 || 2 || 3);
   }

   //Write the output file
   string outputfilename;
   cin >> outputfilename;

   ofstream outputfile;
   outputfile.open(outputfilename.c_str());
   assert(outputfile.is_open());

   outputfile << N_points << endl;
   outputfile << scientific << setprecision(16);
   for(int i_loop=0; i_loop < N_points; i_loop++){
      outputfile << x_coord[i_loop] << "\t" << y_coord[i_loop] << "\t" << z_coord[i_loop] << endl;
   }
   outputfile.close();

   return 0;
}

void calculate_points(int N_points, double radius, vector<double>& coord1, vector<double>& coord2, vector<double>& coord3)
{   
   double small_offset = radius/1e12;

   for(int i_loop=0; i_loop < N_points; i_loop++)
   {
      coord1[i_loop] = radius*cos(2*pi*i_loop/(N_points-1));
      coord2[i_loop] = radius*sin(2*pi*i_loop/(N_points-1));
      coord3[i_loop] = small_offset;
   }
}
