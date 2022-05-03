#include "util.h" //подключили файл заголовка

using namespace std;

void print_vector(const vector<double>& vec) {
	for (auto x : vec) {
		cout << x << endl;
	}
}

double sum_of_el(const vector<double>& vec) {
	double sum = 0;
	for (auto x : vec) {
		sum += x;
	}
	return sum;
}

void read_parameters(Parameters& params) {

	ifstream input("parameters.dat");
	string str;
	vector<double> vec_l; //вектор с параметрами слева
	vector<double> vec_r; //вектор с параметрами справа

	getline(input, str); //считали заголовок
	getline(input, str); //считали заголовок zone 1
	for (int i = 0; i < 7; i++) {
		getline(input, str);
		vec_l.push_back(stod(str.substr(10)));
	}
	getline(input, str); //считали заголовок zone 2
	for (int i = 0; i < 7; i++) {
		getline(input, str);
		vec_r.push_back(stod(str.substr(10)));
	}
	params.P_L = vec_l[0] * params.P0;    //заполнили поля params параметрами слева
	params.T_L = vec_l[1];
	params.V_L = vec_l[2];
	params.C_L_1 = vec_l[3];
	params.C_L_2 = vec_l[4];
	params.C_L_3 = vec_l[5];
	params.C_L_4 = vec_l[6];

	params.P_R = vec_r[0] * params.P0;    //заполнили поля params параметрами справа
	params.T_R = vec_r[1];
	params.V_R = vec_r[2];
	params.C_R_1 = vec_r[3];
	params.C_R_2 = vec_r[4];
	params.C_R_3 = vec_r[5];
	params.C_R_4 = vec_r[6];

	params.parameters_to_vector();
}

void calc_concentration(Parameters& params) {

	double denominator_l = 0, denominator_r = 0; // знаменателm дроби, в котором рассчитываются концентрации. Общий для всех концентраций с одной стороны.

	denominator_l = params.M_C2H2 * params.C_L[0] + params.M_O2 * params.C_L[1] + params.M_N2 * params.C_L[2] + params.M_Ar * params.C_L[3];
	denominator_r = params.M_C2H2 * params.C_R[0] + params.M_O2 * params.C_R[1] + params.M_N2 * params.C_R[2] + params.M_Ar * params.C_R[3];

	for (int i = 0; i < int(params.C_L.size()); i++) {
		params.n_L[i] = params.C_L[i] / denominator_l;
	}
	for (int i = 0; i < int(params.C_R.size()); i++) {
		params.n_R[i] = params.C_R[i] / denominator_r;
	}
}

void calc_adiabatic_index(Parameters& params) {
	double numerator_l = 0; //числители дроби, в которой рассчитывается показатель адиабаты
	double numerator_r = 0;

	numerator_l = params.n_L[0] * params.Cp_H2 + params.n_L[1] * params.Cp_O2 + params.n_L[2] * params.Cp_N2 + params.n_L[3] * params.Cp_Ar;
	numerator_r = params.n_R[0] * params.Cp_H2 + params.n_R[1] * params.Cp_O2 + params.n_R[2] * params.Cp_N2 + params.n_R[3] * params.Cp_Ar;

	params.gl = numerator_l / (numerator_l - params.R0 * sum_of_el(params.n_L));
	params.gr = numerator_r / (numerator_r - params.R0 * sum_of_el(params.n_R));
	cout << params.gl << " , " << params.gr << endl;
}

void calc_start_density(Parameters& params) {
	params.D_L = params.P_L / (params.R0 * params.T_L * sum_of_el(params.n_L));
	params.D_R = params.P_R / (params.R0 * params.T_R * sum_of_el(params.n_R));
}

void calc_sound_velocity(Parameters& params) {
	cout << "params.D_L= " << params.D_L << endl;
	cout << "params.D_R= " << params.D_R << endl;
	params.sv_L = sqrt(params.gl * params.P_L / params.D_L);
	params.sv_R = sqrt(params.gr * params.P_R / params.D_R);
}

void de_dimens(Parameters& params) {
	params.P_L = params.P_L / params.Q;				//обезразмерили начальное давление слева
	params.P_R = params.P_R / params.Q;				//обезразмерили начальное давление справа
	params.V_L = params.V_L / sqrt(params.Q);		//обезразмерили скорость слева
	params.V_R = params.V_R / sqrt(params.Q);		//обезразмерили скорость справа
	params.sv_L = params.sv_L / sqrt(params.Q);		//обезразмерили скорость звука слева
	params.sv_R = params.sv_R / sqrt(params.Q);		//обезразмерили скорость звука слева
}

void re_dimens(Parameters& params, double& s, double& P, double& V, double& E) {
	P = P * params.Q;
	s = s * sqrt(params.Q);
	V = V * sqrt(params.Q);
	E = E * params.Q;
}

double pressure_initial_guess(Parameters& params) {
	double PL = params.P_L;							//переобозначения для удобства записи формулы для начального приближения давления
	double DL = params.D_L;
	double VL = params.V_L;
	double CL = params.sv_L;
	double PR = params.P_R;
	double DR = params.D_R;
	double VR = params.V_R;
	double CR = params.sv_R;
	if ((PL * DR * CR + PR * DL * CL + (VL - VR) * DL * CL * DR * CR) / (DL * CL + DR * CR) < 0) {
		cout << "calc_contact_pressure_velocity --> initial pressure guess is negative" << endl;
	}
	return (PL * DR * CR + PR * DL * CL + (VL - VR) * DL * CL * DR * CR) / (DL * CL + DR * CR);
}

void calc_F_and_DF(double P, double D, double V, double C, double g, double curr_press, double& F, double& DF) {
	double p_ratio;			//отношение давления с предудыщей итерации к актуальному    

	double expr1;			//выражения для удобства записи формул
	double expr2;
	p_ratio = curr_press / P;
	cout << "curr_press= " << curr_press << endl;
	if (curr_press <= P) {
		/* волна разрежения */
		expr1 = 2.0 / (g - 1.0);
		F = expr1 * C * (pow(p_ratio, 1.0 / (expr1 * g)) - 1.0);
		DF = (1.0 / (D * C)) * pow(p_ratio, -0.5 * (g + 1.0) / g);
	}
	else {
		expr2 = sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
		F = (curr_press - P) / (C * D * expr2);
		DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / (g * D * C * pow(expr2, 3.0));
	}
}


void calc_pressure_velocity(Parameters& params, double& P_cont, double& V_cont) {
	double P_old;			//давление на предыдущей итерации
	double fl, fr;			//значения функции слева и справа
	double fld, frd;		//значение производной функции слева и справа
	int iter_num = 0;		// количество проведенных итераций 
	double criteria;		// переменная для определения сходимости 

	/*введение обозначений для удобства записи формул*/
	double PL = params.P_L;
	double DL = params.D_L;
	double VL = params.V_L;
	double CL = params.sv_L;
	double gl = params.gl;
	double PR = params.P_R;
	double DR = params.D_R;
	double VR = params.V_R;
	double CR = params.sv_R;
	double gr = params.gr;

	if ((2.0 * CL / (gr - 1.0)) + (2.0 * CL / (gl - 1.0)) <= VR - VL) {
		/* случай возникновения вакуума */
		cout << "vacuum is generated.rarefaction wave to the left and right." << endl;
	}

	/* расчет начального приближения для давления */
	P_old = pressure_initial_guess(params);
	if (P_old < 0.0) {
		cout << "calc_contact_pressure_velocity --> initial pressure guess is negative" << endl;
		//exit(EXIT_FAILURE);
	}

	/* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
	do {
		calc_F_and_DF(PL, DL, VL, CL, gl, P_old, fl, fld);
		calc_F_and_DF(PR, DR, VR, CR, gr, P_old, fr, frd);
		/*
		cout << "fl= " << fl << endl;
		cout << "fr= " << fr << endl;
		cout << "fld= " << fld << endl;
		cout << "frd= " << frd << endl;
		*/
		P_cont = P_old - (fl + fr + VR - VL) / (fld + frd);
		criteria = 2.0 * fabs((P_cont - P_old) / (P_cont + P_old));
		iter_num++;
		if (iter_num > MAX_ITER_NUM) {
			cout << "calc_contact_pressure_velocity --> number of iterations exceeds the maximum value" << endl;
			exit(EXIT_FAILURE);
		}
		if (P_cont < 0.0) {
			cout << "ncalc_contact_pressure_velocity --> pressure is negative" << endl;
			exit(EXIT_FAILURE);
		}
		P_old = P_cont;
	} while (criteria > params.eps);

	/* скорость контактного разрыва */
	V_cont = 0.5 * (VL + VR + fr - fl);

}

int sample_solid_solution(Parameters& params, double P_cont, double V_cont, double s, double& P_out, double& D_out, double& V_out, double& E_out, bool& flg1, bool& flg2, bool& flg3, bool& flg4) {
	/*введение обозначений для удобства записи формул*/
	double PL = params.P_L;
	double DL = params.D_L;
	double VL = params.V_L;
	double CL = params.sv_L;
	double gl = params.gl;
	double PR = params.P_R;
	double DR = params.D_R;
	double VR = params.V_R;
	double CR = params.sv_R;
	double gr = params.gr;

	double g1l, g2l, g3l, g4l, g5l, g6l, g7l;		//вспомогательные выражения с показателем адиабаты для удобства записи формул
	double g1r, g2r, g3r, g4r, g5r, g6r, g7r;

	/* скорости левых волн */
	double shl, stl;								// скорости "головы" и "хвоста" левой волны разрежения 
	double sl;										// скорость левой ударной волны 

	/* скорости правых волн */
	double shr, str;								// скорости "головы" и "хвоста" правой волны разрежения 
	double sr;										// скорость правой ударной волны 

	double cml, cmr;								// скорости звука слева и справа от контактного разрыва 
	double c;										// локальная скорость звука внутри волны разрежения 
	double p_ratio;
	double d = 0, v = 0, p = 0, e = 0;									// отобранные значения плотности, скорости и давления 

	g1l = 0.5 * (gl - 1.0) / gl;
	g2l = 0.5 * (gl + 1.0) / gl;
	g3l = 2.0 * gl / (gl - 1.0);
	g4l = 2.0 / (gl - 1.0);
	g5l = 2.0 / (gl + 1.0);
	g6l = (gl - 1.0) / (gl + 1.0);
	g7l = 0.5 * (gl - 1.0);

	g1r = 0.5 * (gr - 1.0) / gr;
	g2r = 0.5 * (gr + 1.0) / gr;
	g3r = 2.0 * gr / (gr - 1.0);
	g4r = 2.0 / (gr - 1.0);
	g5r = 2.0 / (gr + 1.0);
	g6r = (gr - 1.0) / (gr + 1.0);
	g7r = 0.5 * (gr - 1.0);


	/*случай возникновения ваккуума*/
	if (VL - VR <= -((2.0 * CR / (gr - 1.0)) + (2.0 * CL / (gl - 1.0)))) {
		/* левая волна разрежения */
		shl = VL - CL;
		cml = 0;
		stl = (2.0 * CL / (gr - 1.0)) + (2.0 * CL / (gl - 1.0));
		shr = VR + CR;
		cmr = 0;
		str = (2.0 * CL / (gr - 1.0)) + (2.0 * CL / (gl - 1.0));
		if (s <= shl) {

			/* параметры слева от разрыва */
			d = DL;
			v = VL;
			p = PL;
			e = (p / (gl - 1)) * params.R0;
			P_out = p;			//выходной результат
			D_out = d;
			V_out = v;
			E_out = e;
			return 1;
		}
		if (s > stl && s < V_cont) {
			/* параметры слева от контактного разрыва */
			d = 0;
			v = 0;
			p = 0;
			e = 0;
			P_out = p;			//выходной результат
			D_out = d;
			V_out = v;
			E_out = e;
			return 1;
		}
		if(s < stl && s >= shl) {
			/* параметры внутри левой волны разрежения */
			v = g5l * (CL + g7l * VL + s);
			c = g5l * (CL + g7l * (VL - s));
			d = DL * pow(c / CL, g4l);
			p = PL * pow(c / CL, g3l);
			e = (p / (gl - 1)) * params.R0;
			if (p > params.eps && d > params.eps)
			{
				P_out = p;			//выходной результат
				D_out = d;
				V_out = v;
				E_out = e;
			}
			else {
				P_out = 0;			//выходной результат
				D_out = 0;
				V_out = v;
				E_out = 0;
			}
			return 1;
		}
		/* правая волна разрежения */
		
		if (s >= shr) {
			/* параметры справа от разрыва */
			d = DR;
			v = VR;
			p = PR;
			e = (p / (gl - 1)) * params.R0;
			P_out = p;			//выходной результат
			D_out = d;
			V_out = v;
			E_out = e;
			return 1;
		}
		if (s <= str && s >= V_cont) {
			/* параметры справа от контактного разрыва */
			d = 0;
			v = 0;
			p = 0;
			e = 0;
			P_out = p;			//выходной результат
			D_out = d;
			V_out = v;
			E_out = e;
			return 1;
		}
		if (s > str && s < shr) {
			/* параметры внутри правой волны разрежения */
			v = g5r * (-CR + g7r * VR + s);
			c = g5r * (CR - g7r * (VR - s));
			d = DR * pow(c / CR, g4r);
			p = PR * pow(c / CR, g3r);
			e = (p / (gl - 1)) * params.R0;
			if (p > params.eps && d > params.eps)
			{
				P_out = p;			//выходной результат
				D_out = d;
				V_out = v;
				E_out = e;
			}
			else {
				P_out = 0;			//выходной результат
				D_out = 0;
				V_out = v;
				E_out = 0;
			}
			return 1;
		}
	}

	if (s <= V_cont) {
		/* рассматриваемая точка - слева от контактного разрыва */
		if (P_cont <= PL) {
			if (flg1 == 0) {
				cout << "rarefaction wave to the left" << endl;
				flg1 = 1;
			}
			/* левая волна разрежения */
			shl = VL - CL;
			if (s <= shl) {
				/* параметры слева от разрыва */
				d = DL;
				v = VL;
				p = PL;
				e = (p / (gl - 1)) * params.R0;
			}
			else {
				cml = CL * pow(P_cont / PL, g1l);
				stl = V_cont - cml;
				if (s > stl) {
					/* параметры слева от контактного разрыва */
					d = DL * pow(P_cont / PL, 1.0 / gl);
					v = V_cont;
					p = P_cont;
					e = (p / (gl - 1)) * params.R0;
				}
				else {
					/* параметры внутри левой волны разрежения */
					v = g5l * (CL + g7l * VL + s);
					c = g5l * (CL + g7l * (VL - s));
					d = DL * pow(c / CL, g4l);
					p = PL * pow(c / CL, g3l);
					e = (p / (gl - 1)) * params.R0;
					if (p > params.eps && d > params.eps)
					{
						P_out = p;			//выходной результат
						D_out = d;
						V_out = v;
						E_out = e;
						return 2;
					}
					else {
						P_out = 0;			//выходной результат
						D_out = 0;
						V_out = v;
						E_out = 0;
						return 2;
					}
				}
			}
		}
		else {
			if (flg2 == 0) {
				cout << "shock wave to the left" << endl;
				flg2 = 1;
			}
			/* левая ударная волна */
			p_ratio = P_cont / PL;
			sl = VL - CL * sqrt(g2l * p_ratio + g1l);
			if (s <= sl) {
				/* параметры слева от разрыва */
				d = DL;
				v = VL;
				p = PL;
				e = (p / (gl - 1)) * params.R0;
			}
			else {
				/* параметры за левой ударной волной */
				d = DL * (p_ratio + g6l) / (p_ratio * g6l + 1.0);
				v = V_cont;
				p = P_cont;
				e = (p / (gl - 1)) * params.R0;
			}
		}
	}
	else {
		if (P_cont > PR) {
			if (flg3 == 0) {
				cout << "shock wave to the right" << endl;
				flg3 = 1;
			}
			/* правая ударная волна */
			p_ratio = P_cont / PR;
			sr = VR + CR * sqrt(g2r * p_ratio + g1r);
			if (s >= sr) {
				/* параметры справа от разрыва */
				d = DR;
				v = VR;
				p = PR;
				e = (p / (gr - 1)) * params.R0;
			}
			else {
				/* параметры за правой ударной волной */
				d = DR * (p_ratio + g6r) / (p_ratio * g6r + 1.0);
				v = V_cont;
				p = P_cont;
				e = (p / (gr - 1)) * params.R0;
			}
		}
		else {
			/* правая волна разрежения */
			if (flg4 == 0) {
				cout << "rarefaction wave to the right" << endl;
				flg4 = 1;
			}
			shr = VR + CR;
			if (s >= shr) {
				/* параметры справа от разрыва */
				d = DR;
				v = VR;
				p = PR;
				e = (p / (gr - 1)) * params.R0;
			}
			else {
				cmr = CR * pow(P_cont / PR, g1r);
				str = V_cont + cmr;
				if (s <= str) {
					/* параметры справа от контактного разрыва */
					d = DR * pow(P_cont / PR, 1.0 / gr);
					v = V_cont;
					p = P_cont;
					e = (p / (gr - 1)) * params.R0;
				}
				else {
					/* параметры внутри правой волны разрежения */
					v = g5r * (-CR + g7r * VR + s);
					c = g5r * (CR - g7r * (VR - s));
					d = DR * pow(c / CR, g4r);
					p = PR * pow(c / CR, g3r);
					e = (p / (gr - 1)) * params.R0;
					if (p > params.eps && d > params.eps)
					{
						P_out = p;			//выходной результат
						D_out = d;
						V_out = v;
						E_out = e;
						return 2;
					}
					else {
						P_out = 0;			//выходной результат
						D_out = 0;
						V_out = v;
						E_out = 0;
						return 2;
					}
				}
			}
		}
	}


	P_out = p;			//выходной результат
	D_out = d;
	V_out = v;
	E_out = e;

	return 3
		;
}

int calc_boundary_velocity(const Parameters& params, const double& P_cont, const double& V_cont, double& BV_L, double& BV_R) {
	double PL = params.P_L;
	double VL = params.V_L;
	double CL = params.sv_L;
	double gl = params.gl;
	double PR = params.P_R;
	double VR = params.V_R;
	double CR = params.sv_R;
	double gr = params.gr;

	double g1l, g2l;
	double g1r, g2r;

	double p_ratio;

	g1l = 0.5 * (gl - 1.0) / gl;
	g2l = 0.5 * (gl + 1.0) / gl;

	g1r = 0.5 * (gr - 1.0) / gr;
	g2r = 0.5 * (gr + 1.0) / gr;


	if (VL - VR < -((2.0 * CR / (gr - 1.0)) + (2.0 * CL / (gl - 1.0)))) {
		/* левая волна разрежения */
		BV_L = VL - CL;

		/* правая волна разрежения */
		BV_R = VR + CR;
		return 1;
	}




	/* рассматриваемая точка - слева от контактного разрыва */
	if (P_cont <= PL) {
		/* левая волна разрежения */
		BV_L = VL - CL;
	}
	else {
		/* левая ударная волна */
		p_ratio = P_cont / PL;
		BV_L = VL - CL * sqrt(g2l * p_ratio + g1l);
	}

	if (P_cont > PR) {
		/* правая ударная волна */
		p_ratio = P_cont / PR;
		BV_R = VR + CR * sqrt(g2r * p_ratio + g1r);

	}
	else {
		/* правая волна разрежения */
		BV_R = VR + CR;
	}
	return 2;
}