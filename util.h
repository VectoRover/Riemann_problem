#pragma once
#include <iostream>	//����������� ������ �����/������
#include <fstream>	//����������� ������ ��� ������ � ������ (��������� iostream)	
#include <string>	//����������� ����� ���������, � ������� �������� ����� string
#include <vector>	//����������� ����� ���������, � ������� �������� ����� vector
#include <cmath>	//����������� ����� ��������� ��� ������ � ��� ��������� (pow, sqrt � ��)
#include "parameters.h" //���������� ���� ���������, � ������� �������� ����� Parameters
#define MAX_ITER_NUM            20

using namespace std;

/*
print_vector �������� � ������� ������ ������.���������� � ����� utils.cpp.
*/
void print_vector(const vector<double>& vec);

/*
sum_of_el ���������� ����� ��������� ������� �� ������. ���������� � ����� utils.cpp.
*/
double sum_of_el(const vector<double>& vec);

/*
read_parameters ��������� ���� parameters.dat, ������ �� ���� ������ � ����������
� ��������������� ���� ������ Parameters. ���������� � ����� utils.cpp.
*/
void read_parameters(Parameters& params);

/*
calc_concentration �� ��������� (���������) ������� ���� ���������� ������������ ���������.
����� ����� �������� ����� ����������� �����, �� ������ � ����������. ���������� � ����� utils.cpp.
*/
void calc_concentration(Parameters& params);

/*
calc_adiabatic_index �� ������������, ������������ ��� ���������� ����������� �
������������� ������� ���������� ������������ ���������� ��������.
���������� � ����� utils.cpp
*/
void calc_adiabatic_index(Parameters& params);

/*
calc_start_density ���������� ���������. ���������� � ����� utils.cpp.
*/
void calc_start_density(Parameters& params);

/*
calc_sound_velocity ���������� �������� ����� ����� � ������ �� ���������� ������� � �������� � ���� sv_L, sv_R ������ Parameters.
���������� � ����� utils.cpp.
*/
void calc_sound_velocity(Parameters& params);

/*
de_dimens ��������������� ��������� ������.
���������� � ����� utils.cpp.
*/
void de_dimens(Parameters& params);

/*
re_dimens ���������� ����������� ���������� ������ ����� �� ������� � ����
*/
void re_dimens(Parameters& params, double& s, double& P, double& V, double& E);

/*
pressure_initiasl_guess ���������� ��������� ����������� �������� ��� ������� �������� �� ��
���������� � ����� utils.cpp.
*/
double pressure_initial_guess(Parameters& params);

/*
calc_F_and_DF ������������ ������� F, ������������ �������� ���� �� ���������� �������, � �� ����������� �� �������� ����� DF
*/
void calc_F_and_DF(double P, double D, double V, double C, double g, double curr_press, double& F, double& DF);

/*
calc_pressure_velocity ������������ �������� � �������� �� ��
*/
void calc_pressure_velocity(Parameters& params, double& P_cont, double& v_cont);


/*
    sample_solid_solution -- ������� ������ �������
    p_cont - �������� �� ���������� �������
    v_cont - �������� �� ���������� �������
    s - �������� x/t, ��� �������� ���������� �������

    v_ncons_res[M] - ������ ���������������� ���������� � ����������� ������������ (out)
*/
int sample_solid_solution(Parameters& params, double P_cont, double V_cont, double s, double& P_out, double& D_out, double& V_out, double& E_out, bool& flg1, bool& flg2, bool& flg3, bool& flg4);


/*
calc_waves_boundary ���������� �������� ��������� ������������������ ���� (�����, ����� ����������
�������� ��������� ������������� ����������)
*/
int calc_boundary_velocity(const Parameters & params, const double& P_cont, const double& V_const, double& BV_L, double& BV_R);
