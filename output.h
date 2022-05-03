#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

ofstream create_header_to_output_file();
void output_to_file(ofstream& output, const double& s, const double& P, const double& D, const double& V, const double& E);
