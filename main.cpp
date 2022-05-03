#include <iostream>
#include <cmath>
#include "util.h"
#include "grid.h"
#include "output.h"
#include "parameters.h" //подключили файл заголовка, в котором хранится класс Parameters

using namespace std;

int main(void) {
	Parameters params;
	double v_cont, p_cont;
	vector<double> xc(CELLS_NUMBER);
	vector<double> x(CELLS_NUMBER + 1);
	double s;	//текущее значение автомодельной переменной
	double P_out = 0, D_out = 0, V_out = 0, E_out = 0; //значения давления, плотности и скорости при отборе решения
	ofstream output("exact_solution.dat");
	double BV_L = 0, BV_R = 0;
	bool flg1 = 0, flg2 = 0, flg3 = 0, flg4 = 0;

	cout << "***   Start programm   ***" << endl;

	read_parameters(params);									//

	//de_dimens(params);		

	calc_concentration(params);									//
	calc_adiabatic_index(params);								//
	calc_start_density(params);
    de_dimens(params);
	calc_sound_velocity(params);								//		

	cout << params.sv_L << " " << params.V_L << endl;
	cout << params.sv_R << " " << params.V_R << endl;
	//params.print_parameters();


	calc_pressure_velocity(params, p_cont, v_cont);
	calc_boundary_velocity(params, p_cont, v_cont, BV_L, BV_R);
	//BV_L = BV_L * sqrt(params.Q);
	//BV_R = BV_R * sqrt(params.Q);
	cout << BV_L << " " << BV_R << endl;
	build_grid(BV_L - abs(BV_L/10.0), BV_R + BV_R/10.0, CELLS_NUMBER, xc, x);


	//output = create_header_to_output_file();

	//цикл по ячейкам
	for (int i_cell = 0; i_cell < CELLS_NUMBER; i_cell++) {
		// изначально решение строится на отрезке [-0.5;0.5] 
		//s = xc[i_cell] / sqrt(params.Q);
		s = xc[i_cell];
		// отбор решения для заданного значения s 
		sample_solid_solution(params, p_cont, v_cont, s, P_out, D_out, V_out, E_out, flg1, flg2, flg3, flg4);
		//возвращение размерностей
		re_dimens(params, s, P_out, V_out, E_out);
		// запись вектора в файл
		output_to_file(output, s * sqrt(params.Q), P_out, D_out, V_out, E_out);
	}

	//params.print_parameters();  
	cout << "***   Finish programm   ***" << endl;


	system("pause");
	return 0;
}
