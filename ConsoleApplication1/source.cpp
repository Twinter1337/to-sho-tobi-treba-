#include <iostream>
#include <cmath>
#include "Header.h"
using namespace std;

void print_matr(double** matr, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			cout << matr[i][j] << "\t";
		}
		cout << endl;
	}
}

double* Find_Coeff(double** pol_def, double** free_member_col_def, int row, int col) {
	double** matr_N = MakeN(pol_def, row, col);

	cout << endl;
	print_matr(matr_N, col, col);
	cout << endl;

	double** N_free_member_col = MakeFreeColN(pol_def, free_member_col_def, row, col);
	cout << endl;
	print_matr(N_free_member_col, col, 1);
	cout << endl;

	double line_det = 0;
	double square_det = 0;
	double cubic_det = 0;

	switch (col) {
	case(2):
		line_det = matr_N[0][0] * matr_N[1][1] - matr_N[0][1] * matr_N[1][0];
		cout << line_det << endl << endl;
		break;
	case(3):
		square_det = DetM(matr_N, col);
		cout << square_det << endl << endl;
		break;
	case(4):
		cubic_det = DetM(matr_N, col);
		cout << cubic_det << endl << endl;
		break;
	}

	double* matr_det_cram = new double[col];
	double* matr_cram_res = new double[col];

	double** cram_matr = MakeMatr(col, col);
	while (true) {
		for (int i = 0; i < col; i++) {
			for (int j = 0; j < col; j++) {
				cram_matr[i][j] = matr_N[i][j];
			}
			cram_matr[i][0] = N_free_member_col[i][0];
		}
		matr_det_cram[0] = DetM(cram_matr, col);
		for (int i = 0; i < col; i++) {
			for (int j = 0; j < col; j++) {
				cram_matr[i][j] = matr_N[i][j];
			}
			cram_matr[i][1] = N_free_member_col[i][0];
		}
		matr_det_cram[1] = DetM(cram_matr, col);

		if (col == 2) { break; }

		for (int i = 0; i < col; i++) {
			for (int j = 0; j < col; j++) {
				cram_matr[i][j] = matr_N[i][j];
			}
			cram_matr[i][2] = N_free_member_col[i][0];
		}
		matr_det_cram[2] = DetM(cram_matr, col);

		if (col == 3) { break; }

		for (int i = 0; i < col; i++) {
			for (int j = 0; j < col; j++) {
				cram_matr[i][j] = matr_N[i][j];
			}
			cram_matr[i][3] = N_free_member_col[i][0];
		}
		matr_det_cram[3] = DetM(cram_matr, col);

		if (col == 4) { break; }
	}

	switch (col) {
	case(2):
		for (int i = 0; i < col; i++) {
			matr_cram_res[i] = matr_det_cram[i] / line_det;
		}
		break;
	case(3):
		for (int i = 0; i < col; i++) {
			matr_cram_res[i] = matr_det_cram[i] / square_det;
		}
		break;
	case(4):
		for (int i = 0; i < col; i++) {
			matr_cram_res[i] = matr_det_cram[i] / cubic_det;
		}
		break;
	}



	return matr_cram_res;
}

long double DetMatr3x3(double** matr) {
	long double det = matr[0][0] * matr[1][1] * matr[2][2] + matr[0][1] * matr[1][2] * matr[2][0] + matr[1][0] * matr[2][1] * matr[0][2]
		- matr[0][2] * matr[1][1] * matr[2][0] - matr[0][1] * matr[1][0] * matr[2][2] - matr[1][2] * matr[2][1] * matr[0][0];
	return det;
}

double DetMatr4x4(double** matr, int ind_i, int ind_j) {
	double det4x4, det;
	double** det_matr = MakeMatr(3, 3);
	int p = 0, q = 0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i != ind_i and j != ind_j) {
				det_matr[p][q] = matr[i][j];

				q++;

				if (q == 3) {
					p++;

					q = 0;
				}
			}
		}
	}
	static int i = 1;

	det = DetMatr3x3(det_matr);
	det4x4 = pow(1, ind_i + ind_j + 2) * matr[ind_i][ind_j] * det;

	return det4x4;
}

double DetM(double** matr, int size) {
	long double res = 0;

	switch (size) {
	case(2):
		res = matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0];
		break;
	case(3):
		res = DetMatr3x3(matr);
		break;
	case(4):
		res = DetMatr4x4(matr, 0, 0) + DetMatr4x4(matr, 0, 1) + DetMatr4x4(matr, 0, 2) + DetMatr4x4(matr, 0, 3);
		break;
	}

	return res;
}

double** MakeMatr(int row, int col) {
	double** matr = new double* [row];
	for (int i = 0; i < row; i++) {
		matr[i] = new double[col];
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			matr[i][j] = 0;
		}
	}

	return matr;
}

double** TransMatr(double** matr, int row, int col) {
	double** TMatr = MakeMatr(col, row);

	for (int i = 0; i < col; i++) {
		for (int j = 0; j < row; j++) {
			TMatr[i][j] = matr[j][i];
		}
	}

	return TMatr;
}

double** MakeN(double** matr, int row, int col) {
	double** TMatr = TransMatr(matr, row, col);
	double** N = MakeMatr(col, col);

	for (int i = 0; i < col; i++) {
		for (int j = 0; j < col; j++) {
			for (int k = 0; k < row; k++) {
				N[i][j] += TMatr[i][k] * matr[k][j];
			}
		}
	}

	return N;
}

double** MakeFreeColN(double** matr, double** res, int row, int col) {
	double** Nres = MakeMatr(col, 1);
	double** TMatr = TransMatr(matr, row, col);

	for (int i = 0; i < col; i++) {
		for (int j = 0; j < col; j++) {
			for (int k = 0; k < row; k++) {
				Nres[i][j] += TMatr[i][k] * res[k][j];
			}
		}
	}

	return Nres;
}

void DelMatr(double** matr, int row) {
	for (int j = 0; j < row; j++) {
		delete matr[j];
	}
	delete matr;
}
