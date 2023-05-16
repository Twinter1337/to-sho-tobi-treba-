#include <iostream>
#include <cmath>
#include "Header.h"
using namespace std;

int main(void) {
	double line_pol_def[6][2] = { {1.0, 0.41},
								  {1.0, 0.46},
								  {1.0, 0.52},
								  {1.0, 0.61},
								  {1.0, 0.66},
								  {1.0, 0.73} };
	double square_pol_def[6][3] = { {1.0, 0.41,pow(0.41,2)},
									{1.0, 0.46,pow(0.46,2)},
									{1.0, 0.52,pow(0.52,2)},
									{1.0, 0.61,pow(0.61,2)},
									{1.0, 0.66,pow(0.66,2)},
									{1.0, 0.73,pow(0.73,2)} };
	double cubic_pol_def[6][4] = { {1.0, 0.41,pow(0.41,2),pow(0.41,3)},
									{1.0, 0.46,pow(0.46,2),pow(0.46,3)},
									{1.0, 0.52,pow(0.52,2),pow(0.52,3)},
									{1.0, 0.61,pow(0.61,2),pow(0.61,3)},
									{1.0, 0.66,pow(0.66,2),pow(0.66,3)},
									{1.0, 0.73,pow(0.73,2),pow(0.73,3)} };
	double free_member_col_def[6][1] = { {2.73},
										 {2.31},
										 {1.97},
										 {1.76},
										 {1.53},
										 {1.31} };
	double** line_pol = MakeMatr(6, 2);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 2; j++) {
			line_pol[i][j] = line_pol_def[i][j];
		}
	}
	double** square_pol = MakeMatr(6, 3);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 3; j++) {
			square_pol[i][j] = square_pol_def[i][j];
		}
	}
	double** cubic_pol = MakeMatr(6, 4);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++) {
			cubic_pol[i][j] = cubic_pol_def[i][j];
		}
	}
	double** free_mem_col = MakeMatr(6, 1);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 1; j++) {
			free_mem_col[i][j] = free_member_col_def[i][j];
		}
	}

	double* line = Find_Coeff(line_pol, free_mem_col, 6, 2);
	for (int i = 0; i < 2; i++) {
		cout << line[i] << " ";
	}

	DelMatr(line_pol, 6);

	cout << endl;
	double* square = Find_Coeff(square_pol, free_mem_col, 6, 3);
	for (int i = 0; i < 3; i++) {
		cout << square[i] << " ";
	}

	DelMatr(square_pol, 6);

	cout << endl;
	double* cubic = Find_Coeff(cubic_pol, free_mem_col, 6, 4);
	for (int i = 0; i < 4; i++) {
		cout << cubic[i] << " ";
	}

	return 0;
}