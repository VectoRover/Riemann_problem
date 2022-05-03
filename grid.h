#pragma once
#include <vector>
#include "parameters.h"

#define CELLS_NUMBER 1000
#define LEFT_BOUNDARY -0.2 
#define RIGHT_BOUNDARY 0.2
#define STOP_TIME 0.00001

using namespace std;

/* ����������� ��������� ����� � ������� ����� �����

   left_boundary_x - ���������� ����� ������� ��������� ������� (in)
   right_boundary_x - ���������� ������ ������� ��������� ������� (in)
   cells_num - ���������� ����� ����� (in)

   *xc - ������ ��������� ������� ����� (out)
   *x - ������ ��������� ����� ����� (out) */
void build_grid(double left_boundary_x, double right_boundary_x, int cells_num, vector<double> & xc, vector<double> & x); 
