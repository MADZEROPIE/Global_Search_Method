
#ifndef _MAP
#define _MAP

#include <stdlib.h>

//����������� ����� x �� ����������� ������� [-0.5, 0.5] � ��������� 2^(-m), key = 1
void mapd (double x, int m, double* y, int n, int key);    /* map x to y         */

//��� ��� ������� ���� ������������ �� ����
void invmad ( int, double *, int , int *, double *, int , int );  /* map y to x         */
void xyd ( double *, int, double *, int );        /* get preimage       */

#endif