#ifndef _CLIPPER_REGION_H_
#define _CLIPPER_REGION_H_
#pragma once
#include <math.h>
#include <algorithm>
//#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
//#include <iomanip>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <string>

#include "clipper_.h"
#include "clipper_offset.h"
//#include "svg.h"


using namespace std;
using namespace clipperlib;
//using namespace svglib;

double clipper_intersection_area(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double upperleft_x, double upperleft_y,
	double downright_x, double downright_y);
	
double clipper_intersection_area_D(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y,
	double upperleft_x, double upperleft_y,
	double downright_x, double downright_y);

double clipper_intersection_area_H6(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y,
	double v5_x, double v5_y, double v6_x, double v6_y,
	double upperleft_x, double upperleft_y,
	double downright_x, double downright_y);

double clipper_intersection_area_H5(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y,
	double v5_x, double v5_y,
	double upperleft_x, double upperleft_y,
	double downright_x, double downright_y);
#endif _CLIPPER_REGION_H_




