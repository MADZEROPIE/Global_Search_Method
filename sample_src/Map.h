
#ifndef _MAP
#define _MAP

#include <stdlib.h>

//отображение точки x на многомерную область [-0.5, 0.5] с точностью 2^(-m), key = 1
void mapd (double x, int m, double* y, int n, int key);    /* map x to y         */

//Эти две функции пока использовать не надо
void invmad ( int, double *, int , int *, double *, int , int );  /* map y to x         */
void xyd ( double *, int, double *, int );        /* get preimage       */

#endif