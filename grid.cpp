#include "grid.h"

void build_grid(double left_boundary_x, double right_boundary_x, int cells_num, vector<double>& xc, vector<double>& x) {
    int i;
    double h;   /* шаг сетки */

   /* пока реализовано только равномерное распределение узлов сетки */
    h = (right_boundary_x - left_boundary_x) / cells_num;

    /* координаты узлов */
    for (i = 0; i < cells_num + 1; i++) {
        x[i] = left_boundary_x + i * h;
    }

    /* координаты центров ячеек */
    for (i = 0; i < cells_num; i++) {
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }

}