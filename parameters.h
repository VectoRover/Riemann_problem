#pragma once
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Parameters {

public:
	Parameters() {											//����������� ������ �� ��������� �������� ��� ������������ ���������
		P_L = V_L = T_L = P_R = V_R = T_R 
		= C_L_1 = C_L_2 = C_L_3 = C_L_4 = C_L_5 = C_L_6 = C_L_7 = C_L_8 = C_L_9 = C_L_10 = C_L_11 = C_L_12 = C_L_13 = C_L_14 
		= C_R_1 = C_R_2 = C_R_3 = C_R_4 = C_R_5 = C_R_6 = C_R_7 = C_R_8 = C_R_9 = C_R_10 = C_R_11 = C_R_12 = C_R_13 = C_R_14
		= n_L_1 = n_L_2 = n_L_3 = n_L_4 = n_L_5 = n_L_6 = n_L_7 = n_L_8 = n_L_9 = n_L_10 = n_L_11 = n_L_12 = n_L_13 = n_L_14 
		= n_R_1 = n_R_2 = n_R_3 = n_R_4 = n_R_5 = n_R_6 = n_R_7 = n_R_8 = n_R_9 = n_R_10 = n_R_11 = n_R_12 = n_R_13 = n_R_14
		= gl = gr = D_L = D_R =	sv_L = sv_R = 0;
	}

	double R0 = 8.31;		//������������� ������� ����������,  �� * ����^-1 * K^-1
	
	double M_C2H2 = 0.0020158;		//  �������� �����, ��/����
	double M_O2 = 0.0319988;
	double M_N2 = 0.0280134;
	double M_Ar = 0.039948;
	double M_CO2 = 0.0440098;
	double M_H2O = 0.0180152;
	double M_H = 0.0010079;
	double M_O = 0.0159994;
	double M_OH = 0.0170073;
	double M_HO2 = 0.0330067;
	double M_CO = 0.0280104;
	double M_H2CO = 0.0300262;
	double M_HCO = 0.0290183;
	double M_CH3 = 0.0150347;

	double Cp_C2H2 = 69.84062706;       // �� * ����^-1 * K^-1
	double Cp_O2 = 35.53481848;
	double Cp_N2 = 33.74732673;
	double Cp_Ar = 20.7860066;
	double Cp_CO2 = 53.84930693;
	double Cp_H2O = 45.77349835;
	double Cp_H = 20.78765677;
	double Cp_O = 21.11405941;
	double Cp_OH = 32.83188119;
	double Cp_HO2 = 49.49036304;
	double Cp_CO = 34.0169967;
	double Cp_H2CO = 64.0830033;
	double Cp_HCO = 49.00471947;
	double Cp_CH3 = 62.67224422;

	double Q = 41'868'000.0;		//�������� � ������������ m^2 / s^2 ��� ����������������
	double eps = 1 * pow(10, -7);	//������������ �������� �� �������
	double L = 0.0006470548662;		//�������� � ����������� m ��� ���������������� �������� ����������

	double P0 = 101325;				//��� ��������, ��

	double P_L, V_L, T_L; //��������, �������� � ����������� ����� �� �������. ��������� �� �����
	double P_R, V_R, T_R; //��������, �������� � ����������� ������ �� �������. ��������� �� �����

	double C_L_1, C_L_2, C_L_3, C_L_4, C_L_5, C_L_6, C_L_7, C_L_8, C_L_9, C_L_10, C_L_11, C_L_12, C_L_13, C_L_14; //�������� ��������� ��������� �����. ��������� �� �����
	double C_R_1, C_R_2, C_R_3, C_R_4, C_R_5, C_R_6, C_R_7, C_R_8, C_R_9, C_R_10, C_R_11, C_R_12, C_R_13, C_R_14; //�������� ��������� ��������� ������. ��������� �� �����

	double n_L_1, n_L_2, n_L_3, n_L_4, n_L_5, n_L_6, n_L_7, n_L_8, n_L_9, n_L_10, n_L_11, n_L_12, n_L_13, n_L_14; //������������ ��������� �����. ������������ � ������� calc_concentration
	double n_R_1, n_R_2, n_R_3, n_R_4, n_R_5, n_R_6, n_R_7, n_R_8, n_R_9, n_R_10, n_R_11, n_R_12, n_R_13, n_R_14; //������������ ��������� ������. ������������ � ������� calc_concentration

	double gl; //���������� �������� ����� �����
	double gr; //���������� �������� ����� ������

	double D_L; //���������. ������������ � ������� calc_density
	double D_R; //���������. ������������ � ������� calc_density

	double sv_L; //�������� ����� ����� �� ���������� �������
	double sv_R; //�������� ����� ������ �� ���������� �������

	vector<double> C_L; //�� �� ���������, ��� � ����, �� ���������� ������ ��� �������� �������� � �������
	vector<double> C_R; //�� �� ���������, ��� � ����, �� ���������� ������ ��� �������� �������� � �������
	vector<double> n_L;	//�� �� ���������, ��� � ����, �� ���������� ������ ��� �������� �������� � �������
	vector<double> n_R;	//�� �� ���������, ��� � ����, �� ���������� ������ ��� �������� �������� � �������

	void parameters_to_vector();
	void print_parameters();
private:

};


