#pragma once
#include <iostream>	//подключение потока ввода/вывода
#include <fstream>	//подключение класса для чтения и записи (наследник iostream)	
#include <string>	//подключение файла заголовка, в котором хранится класс string
#include <vector>	//подключение файла заголовка, в котором хранится класс vector
#include <cmath>	//подключение файла заголовка для работы с мат функциями (pow, sqrt и тд)
#include "parameters.h" //подключили файл заголовка, в котором хранится класс Parameters
#define MAX_ITER_NUM            20

using namespace std;

/*
print_vector печатает в консоль вектор даблов.Реализация в файле utils.cpp.
*/
void print_vector(const vector<double>& vec);

/*
sum_of_el возвращает сумму элементов вектора из даблов. Реализация в файле utils.cpp.
*/
double sum_of_el(const vector<double>& vec);

/*
read_parameters открывает файл parameters.dat, читает из него данные и записывает
в соответствующие поля класса Parameters. Реализация в файле utils.cpp.
*/
void read_parameters(Parameters& params);

/*
calc_concentration по молярному (объемному) составу газа определяет концентрации компонент.
Нужно знать молярные массы компонентов смеси, их задаем в реализации. Реализация в файле utils.cpp.
*/
void calc_concentration(Parameters& params);

/*
calc_adiabatic_index по концентрации, теплоемкости при постоянной температуре и
универсальной газовой постоянной рассчитывает показатель адиабаты.
Реализация в файле utils.cpp
*/
void calc_adiabatic_index(Parameters& params);

/*
calc_start_density возвращает плотность. Реализация в файле utils.cpp.
*/
void calc_start_density(Parameters& params);

/*
calc_sound_velocity определяет скорости звука слева и справа от начального разрыва и помещает в поля sv_L, sv_R класса Parameters.
Реализация в файле utils.cpp.
*/
void calc_sound_velocity(Parameters& params);

/*
de_dimens обезразмеривает параметры задачи.
Реализация в файле utils.cpp.
*/
void de_dimens(Parameters& params);

/*
re_dimens возвращает размерности параметрам задачи перед их печатью в файл
*/
void re_dimens(Parameters& params, double& s, double& P, double& V, double& E);

/*
pressure_initiasl_guess возвращает начальное приближение давления для расчета давления на КР
Реализация в файле utils.cpp.
*/
double pressure_initial_guess(Parameters& params);

/*
calc_F_and_DF рассчитывает функцию F, определяющую скорость газа на контактном разрыве, и ее производную по давлению среды DF
*/
void calc_F_and_DF(double P, double D, double V, double C, double g, double curr_press, double& F, double& DF);

/*
calc_pressure_velocity рассчитывает давление и скорость на КР
*/
void calc_pressure_velocity(Parameters& params, double& P_cont, double& v_cont);


/*
    sample_solid_solution -- функция отбора решения
    p_cont - давление на контактном разрыве
    v_cont - скорость на контактном разрыве
    s - значение x/t, для которого отбирается решение

    v_ncons_res[M] - вектор неконсервативных переменных с отобранными компонентами (out)
*/
int sample_solid_solution(Parameters& params, double P_cont, double V_cont, double s, double& P_out, double& D_out, double& V_out, double& E_out, bool& flg1, bool& flg2, bool& flg3, bool& flg4);


/*
calc_waves_boundary возвращает значения скоростей распространяющихся волн (нужно, чтобы определить
диапазон изменения автомодельной переменной)
*/
int calc_boundary_velocity(const Parameters & params, const double& P_cont, const double& V_const, double& BV_L, double& BV_R);
