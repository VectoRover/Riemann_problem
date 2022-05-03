#include "output.h"

using namespace std;

ofstream create_header_to_output_file() {
	ofstream output("exact_solution.dat");
	output << setw(20) << right << "x/t";
	output << setw(20) << right << "P";
	output << setw(20) << right << "D";
	output << setw(20) << right << "V";
	output << setw(20) << right << "E";
	return output;
}

void output_to_file(ofstream& output, const double& s, const double& P, const double& D, const double& V, const double& E) {
	output << setw(20) << right << fixed << setprecision(8) << s;
	output << setw(20) << right << fixed << setprecision(8) << P;
	output << setw(20) << right << fixed << setprecision(8) << D;
	output << setw(20) << right << fixed << setprecision(8) << V;
	output << setw(20) << right << fixed << setprecision(8) << E;
	output << endl;


}