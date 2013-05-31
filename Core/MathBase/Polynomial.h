#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include "defines.h"
#include <stdarg.h>

typedef struct Polynomial {
  ulong m_degree;
  double *m_data;
}Polynomial;

Polynomial * PolynomialInit(const uint degree);
Polynomial * PolynomialCopy(const Polynomial* pol);
Polynomial * PolynomialDelete(Polynomial* pol);

void PolynomialShow(const Polynomial* pol);
void PolynomialSetConstants(Polynomial* pol, ...);
void PolynomialSetConstantsAndDegree(Polynomial* pol, uint degree, ...);
double PolynomialFx(const Polynomial* pol, const double x);
#endif