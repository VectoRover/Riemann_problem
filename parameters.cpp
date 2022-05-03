#include "parameters.h"
#include "util.h"
using namespace std;

void Parameters::print_parameters() {
	cout << "Zone L" << endl;
	cout << P_L << endl;
	cout << V_L << endl;
	cout << T_L << endl;
	cout << C_L_1 << ", " << C_L_2 << ", " << C_L_3 << endl;
	cout << gl << endl;
	cout << D_L << endl;
	cout << sv_L << endl;
	//	print_vector(n_R);

	cout << "*************" << endl;
	cout << "Zone R" << endl;
	cout << P_R << endl;
	cout << V_R << endl;
	cout << T_R << endl;
	cout << C_R_1 << ", " << C_R_2 << ", " << C_R_3 << endl;
	cout << gr << endl;
	cout << D_R << endl;
	cout << sv_R << endl;
	//	print_vector(n_L);

}

void Parameters::parameters_to_vector() {
	C_L.push_back(C_L_1);
	C_L.push_back(C_L_2);
	C_L.push_back(C_L_3);
	C_L.push_back(C_L_4);
	C_L.push_back(C_L_5);
	C_L.push_back(C_L_6);
	C_L.push_back(C_L_7);
	C_L.push_back(C_L_8);
	C_L.push_back(C_L_9);
	C_L.push_back(C_L_10);
	C_L.push_back(C_L_11);
	C_L.push_back(C_L_12);
	C_L.push_back(C_L_13);
	C_L.push_back(C_L_14);

	C_R.push_back(C_R_1);
	C_R.push_back(C_R_2);
	C_R.push_back(C_R_3);
	C_R.push_back(C_R_4);
	C_R.push_back(C_R_5);
	C_R.push_back(C_R_6);
	C_R.push_back(C_R_7);
	C_R.push_back(C_R_8);
	C_R.push_back(C_R_9);
	C_R.push_back(C_R_10);
	C_R.push_back(C_R_11);
	C_R.push_back(C_R_12);
	C_R.push_back(C_R_13);
	C_R.push_back(C_R_14);

	n_L.push_back(n_L_1);
	n_L.push_back(n_L_2);
	n_L.push_back(n_L_3);
	n_L.push_back(n_L_4);
	n_L.push_back(n_L_5);
	n_L.push_back(n_L_6);
	n_L.push_back(n_L_7);
	n_L.push_back(n_L_8);
	n_L.push_back(n_L_9);
	n_L.push_back(n_L_10);
	n_L.push_back(n_L_11);
	n_L.push_back(n_L_12);
	n_L.push_back(n_L_13);
	n_L.push_back(n_L_14);

	n_R.push_back(n_R_1);
	n_R.push_back(n_R_2);
	n_R.push_back(n_R_3);
	n_R.push_back(n_R_4);
	n_R.push_back(n_R_5);
	n_R.push_back(n_R_6);
	n_R.push_back(n_R_7);
	n_R.push_back(n_R_8);
	n_R.push_back(n_R_9);
	n_R.push_back(n_R_10);
	n_R.push_back(n_R_11);
	n_R.push_back(n_R_12);
	n_R.push_back(n_R_13);
	n_R.push_back(n_R_14);
}