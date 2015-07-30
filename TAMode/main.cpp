#include <iostream>
#include <algorithm>
#include <exception>
#include <string>
#include <fstream>
#include <stdexcept>

#include "reactCode.h"
#include "sundials_nvector.h"
#include "nvector_serial.h"
#include "TAMode.h"

#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/accumulators.hpp>

using namespace boost::accumulators;

double times[10] = {0, 20, 40, 60, 80, 100, 120, 160, 200, 240}; 
double Gass = 1.25; //Gas6 dose (in nM)
int kTPS = 10; //number of time points - length of the times array

//lengths need to change with the times[] - 3*length of times[]
double pY[30];  
double tot[30];
double surf[30];

 int gen_num = 71592;  //Number of generations of DREAM sequence
 int par_num = 8;  //Number of parameters fitted via DREAM sequence
 double x[71592][9]; //x[gen_num][par_num+2]
 double Dparams[8];

int main()
{
	//sets up accumulator set for each time point in total
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct1;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct2;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct3;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct4;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct5;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct6;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct7;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct8;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct9;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > acct10;
    
    //sets up accumulator set for each time point in surf
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs1;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs2;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs3;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs4;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs5;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs6;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs7;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs8;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs9;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accs10;
    
    //sets up accumulator set for each time point in pY
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp1;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp2;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp3;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp4;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp5;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp6;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp7;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp8;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp9;
    accumulator_set<double, stats<tag::mean, tag::median, tag::variance> > accp10;

//Deletes previous copies of the output files
  remove("AXL_total_a1_k3.txt");
  remove("AXL_surf_a1_k3.txt");
  remove("AXL_pY_a1_k3.txt");
  
  remove("AXL_total_a1_k3_stats.txt");
  remove("AXL_surf_a1_k3_stats.txt");
  remove("AXL_pY_a1_k3_stats.txt");

   
  std::ifstream in("a1.txt");
  
    //takes values from DREAM text file into holding array
   for (int h = 0; h<gen_num; h++){ 
   	
      for (int i = 0; i<(par_num+1); i++){
   	
	   in>>x[h][i];
   }
}
   
  std::ofstream TAM;
 
   for (int h = 0; h<gen_num; h++){
   
	for (int i = 1; i<(par_num+1); i++){
		//Throws out the error and generation numbers from the holding array, taking one set of parameters
		//at a time into the parameter array
		
		Dparams[i-1] = pow(10, x[h][i]);	
		
	}

	//Dparams[2] = 0;  //Sets AXL expression to 0
	calcProfileMat(Dparams, times, kTPS, tot, pY, surf, Gass);
	//runs the simulation for the current set of parameters
	
	int receptor = 0;
	//Select which receptor data is outputted (0=AXL, 1=Mer, 2=Tyro)
	
	//Adds data points from total, surf, and pY to accumulator sets	
	acct1(tot[receptor + 0]);
	acct2(tot[receptor + 3]);
	acct3(tot[receptor + 6]);
	acct4(tot[receptor + 9]);
	acct5(tot[receptor + 12]);
	acct6(tot[receptor + 15]);
	acct7(tot[receptor + 18]);
	acct8(tot[receptor + 21]);
	acct9(tot[receptor + 24]);
	acct10(tot[receptor + 27]);
	
	accs1(surf[receptor + 0]);
	accs2(surf[receptor + 3]);
	accs3(surf[receptor + 6]);
	accs4(surf[receptor + 9]);
	accs5(surf[receptor + 12]);
	accs6(surf[receptor + 15]);
	accs7(surf[receptor + 18]);
	accs8(surf[receptor + 21]);
	accs9(surf[receptor + 24]);
	accs10(surf[receptor + 27]);
	
	accp1(pY[receptor + 0]);
	accp2(pY[receptor + 3]);
	accp3(pY[receptor + 6]);
	accp4(pY[receptor + 9]);
	accp5(pY[receptor + 12]);
	accp6(pY[receptor + 15]);
	accp7(pY[receptor + 18]);
	accp8(pY[receptor + 21]);
	accp9(pY[receptor + 24]);
	accp10(pY[receptor + 27]);
	
	//Creates text output files for full data sets
	for(int j = 0; j<(kTPS*3); j = j+3){
	//outputs the total at each time point for each parameter set
	TAM.open("AXL_total_a1_k3.txt", std::ofstream::app);
  	TAM<<tot[j+receptor]<<"   ";
  	TAM.close();
  	
  	
  	//outputs the surf at each time point for each parameter set
	TAM.open("AXL_surf_a1_k3.txt", std::ofstream::app);
  	TAM<<surf[j+receptor]<<"   ";
  	TAM.close();
  	
  	//outputs the surf at each time point for each parameter set
	TAM.open("AXL_pY_a1_k3.txt", std::ofstream::app);
  	TAM<<pY[j+receptor]<<"   ";
  	TAM.close();
  	  	
  	}
	//Ends lines to format text file
	TAM.open("AXL_total_a1_k3.txt", std::ofstream::app);
  	TAM<<std::endl;
  	TAM.close();
  	
  	//Ends lines to format text file
	TAM.open("AXL_surf_a1_k3.txt", std::ofstream::app);
  	TAM<<std::endl;
  	TAM.close();
  	
  	//Ends lines to format text file
	TAM.open("AXL_pY_a1_k3.txt", std::ofstream::app);
  	TAM<<std::endl;
  	TAM.close();
  	
  		
}
	
	//Creates and writes to files that give summary statistcs from data sets
	TAM.open("AXL_total_a1_k3_stats.txt", std::ofstream::app);
	TAM<<"Mean: "<<mean(acct1)<<"   "<<mean(acct2)<<"   "<<mean(acct3)<<"   "<<mean(acct4)<<"   "<<mean(acct5)<<"   "<<mean(acct6)<<"   "<<mean(acct7)<<"   "<<mean(acct8)<<"   "<<mean(acct9)<<"   "<<mean(acct10)<<std::endl;
	TAM<<"Median: "<<median(acct1)<<"   "<<median(acct2)<<"   "<<median(acct3)<<"   "<<median(acct4)<<"   "<<median(acct5)<<"   "<<median(acct6)<<"   "<<median(acct7)<<"   "<<median(acct8)<<"   "<<median(acct9)<<"   "<<median(acct10)<<std::endl;
	TAM<<"StdDev: "<<sqrt(variance(acct1))<<"   "<<sqrt(variance(acct2))<<"   "<<sqrt(variance(acct3))<<"   "<<sqrt(variance(acct4))<<"   "<<sqrt(variance(acct5))<<"   "<<sqrt(variance(acct6))<<"   "<<sqrt(variance(acct7))<<"   "<<sqrt(variance(acct8))<<"   "<<sqrt(variance(acct9))<<"   "<<sqrt(variance(acct10))<<std::endl;
	TAM.close();
	
	TAM.open("AXL_surf_a1_k3_stats.txt", std::ofstream::app);
	TAM<<"Mean: "<<mean(accs1)<<"   "<<mean(accs2)<<"   "<<mean(accs3)<<"   "<<mean(accs4)<<"   "<<mean(accs5)<<"   "<<mean(accs6)<<"   "<<mean(accs7)<<"   "<<mean(accs8)<<"   "<<mean(accs9)<<"   "<<mean(accs10)<<std::endl;
	TAM<<"Median: "<<median(accs1)<<"   "<<median(accs2)<<"   "<<median(accs3)<<"   "<<median(accs4)<<"   "<<median(accs5)<<"   "<<median(accs6)<<"   "<<median(accs7)<<"   "<<median(accs8)<<"   "<<median(accs9)<<"   "<<median(accs10)<<std::endl;
	TAM<<"StdDev: "<<sqrt(variance(accs1))<<"   "<<sqrt(variance(accs2))<<"   "<<sqrt(variance(accs3))<<"   "<<sqrt(variance(accs4))<<"   "<<sqrt(variance(accs5))<<"   "<<sqrt(variance(accs6))<<"   "<<sqrt(variance(accs7))<<"   "<<sqrt(variance(accs8))<<"   "<<sqrt(variance(accs9))<<"   "<<sqrt(variance(accs10))<<std::endl;
	TAM.close();
	
	TAM.open("AXL_pY_a1_k3_stats.txt", std::ofstream::app);
	TAM<<"Mean: "<<mean(accp1)<<"   "<<mean(accp2)<<"   "<<mean(accp3)<<"   "<<mean(accp4)<<"   "<<mean(accp5)<<"   "<<mean(accp6)<<"   "<<mean(accp7)<<"   "<<mean(accp8)<<"   "<<mean(accp9)<<"   "<<mean(accp10)<<std::endl;
	TAM<<"Median: "<<median(accp1)<<"   "<<median(accp2)<<"   "<<median(accp3)<<"   "<<median(accp4)<<"   "<<median(accp5)<<"   "<<median(accp6)<<"   "<<median(accp7)<<"   "<<median(accp8)<<"   "<<median(accp9)<<"   "<<median(accp10)<<std::endl;
	TAM<<"StdDev: "<<sqrt(variance(accp1))<<"   "<<sqrt(variance(accp2))<<"   "<<sqrt(variance(accp3))<<"   "<<sqrt(variance(accp4))<<"   "<<sqrt(variance(accp5))<<"   "<<sqrt(variance(accp6))<<"   "<<sqrt(variance(accp7))<<"   "<<sqrt(variance(accp8))<<"   "<<sqrt(variance(accp9))<<"   "<<sqrt(variance(accp10))<<std::endl;
	TAM.close();
	
    return 0;
}
