/**
 * Copyright (c) Quan Bui <mquan.bui@gmail.com>. 
 * An 2D implementation of "Fast Poisson Disk Sampling in Arbitrary Dimensions, R. Bridson, ACM SIGGRAPH 2007 Sketches Program".
 * Anybody could use this code freely, if you feel it good or helpful, please tell me, thank you very much.
 */

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
 
#include <cmath>
#include <ctime>
//#include <chrono>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include "fpd.h"


using namespace std;

int main(int argc, char* argv[])
{
	//Long64 MaxPTS = 10000;
	char* type = new char[128];
	double domain_size, min_dist, max_dist, mean, std_dev;
	unsigned seed;
	int k;
	char* output_name = new char[128];
	char* default_output_name = new char[128];
	strcpy(default_output_name,"tmp_pd_pts.inp");
	domain_size = min_dist = max_dist = 0.0f;

	if (argc < 4) {
		printf("Usage: <FastPDSampling> <seed> <domain_size> <type> <type_arg1> ... <type_argn> <density_param> <output_file> \n");
		printf("Type: <Uniform> <min_dist> \n");
		printf("      <LogNormal> <mean> <std_dev> <max_dist> <min_dist> \n");
		exit(1);
	}
	else {
		//for (int i = 0; i < argc; i++) {
		//	printf("Argument no.%d: %s \n", i, argv[i]);
		//}
		seed = atoi(argv[1]);
		domain_size = atof(argv[2]);
		if (std::string(argv[3]) == "Uniform" && argc > 4) {
			strcpy(type, "Uniform");
			min_dist = atof(argv[4]);
			k = atoi(argv[5]);
			if (argc == 6) {
				strcpy(output_name, default_output_name);
			}
			else {
				strcpy(output_name, argv[6]);
			}			
		}
		else if (std::string(argv[3]) == "LogNormal" && argc > 7) {
			max_dist = atof(argv[4]);
			min_dist = atof(argv[5]);
			mean = atof(argv[6]);
			std_dev = atof(argv[7]);
			k = atoi(argv[8]);
			if (argc == 9) {
				strcpy(output_name, default_output_name);
			}
			else {
				strcpy(output_name, argv[9]);
			}	
		}
		else {
			printf("ERROR!!! Cannot recognize command line argument(s). Please check your command line input. \n");
			printf("Usage: <FastPDSampling> <seed> <domain_size> <type> <type_arg1> ... <type_argn> <density_param> <output_file> \n");
			printf("Type: <Uniform> <min_dist> \n");
			printf("      <LogNormal> <max_dist> <min_dist> <mean> <std_dev> \n");
			exit(1);
		}
	}
	
	srand(seed);
	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//unsigned seed = rand();
	printf("Seed value: %d\n", seed);
	
	if (std::string(argv[3]) == "Uniform") {
		if (domain_size == 0 || min_dist == 0) {
			printf("Error in input values \n");
			exit(1);
		}
		printf("Uniform Distribution \n");
		PDSUniform* ptr_point_set = new PDSUniform(&domain_size, type, &min_dist, &k);
		ptr_point_set->PrintInfo();
		//ptr_point_set->SetSamplePoints(sample_points);
		ptr_point_set->FastPoissonDiskSampling2D();
		printf("Output is written in: %s \n", output_name);
		ptr_point_set->PrintOutput(output_name);
		delete ptr_point_set;
	}
	else if (std::string(argv[3]) == "LogNormal") {
		if (domain_size == 0 || max_dist == 0 || min_dist == 0) {
			printf("Error in input values \n");
			exit(1);
		}
		printf("Log Normal Distribution \n");
		PDSLNorm* ptr_point_set = new PDSLNorm(&domain_size, type, &seed, &max_dist, &min_dist, &mean, &std_dev, &k);
		ptr_point_set->PrintInfo();
		//ptr_point_set->SetSamplePoints(sample_points);
		ptr_point_set->FastPoissonDiskSampling2D();
		printf("Output is written in: %s \n", output_name);
		ptr_point_set->PrintOutput(output_name);
		delete ptr_point_set;
	}
	else {
		printf("ERROR!!! No point set is created, check input and workflow \n");
	}

	return 0;
}

