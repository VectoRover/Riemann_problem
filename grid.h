#pragma once
#include <vector>
#include "parameters.h"

#define CELLS_NUMBER 1000
#define LEFT_BOUNDARY -0.2 
#define RIGHT_BOUNDARY 0.2
#define STOP_TIME 0.00001

using namespace std;

/* Определение координат узлов и центров ячеек сетки

   left_boundary_x - координата левой границы расчетной области (in)
   right_boundary_x - координата правой границы расчетной области (in)
   cells_num - количество ячеек сетки (in)

   *xc - массив координат центров ячеек (out)
   *x - массив координат узлов сетки (out) */
void build_grid(double left_boundary_x, double right_boundary_x, int cells_num, vector<double> & xc, vector<double> & x); 
