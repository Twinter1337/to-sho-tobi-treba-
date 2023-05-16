#ifndef LAB6_NM_HEADER_H
#define LAB6_NM_HEADER_H

void print_matr(double** matr, int row, int col);

double* Find_Coeff(double** pol_matr, double** free_member_col, int row, int col);

long double DetMatr3x3(double** matr);

double DetMatr4x4(double** matr, int ind_i, int ind_j);

double DetM(double** M, int n);

double** MakeMatr(int row, int col);

double** TransMatr(double** matr, int row, int col);

double** MakeN(double** matr, int Size1, int Size2);

double** MakeFreeColN(double** matr, double** res, int row, int col);

void DelMatr(double** matr, int row);

#endif

