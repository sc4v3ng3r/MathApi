#include "Polynomial.h"
#include <math.h>

Polynomial* PolynomialInit(const uint degree)
{
  Polynomial *pol = (Polynomial *) malloc(sizeof(Polynomial));
  if (!pol)
    return NULL;
  
  pol->m_data = (double*) malloc((degree+1) * sizeof(double));
  if (!pol->m_data)
    return NULL; // aqui tem de ter um assert, tratamento de falhas...
  
  pol->m_degree = degree;
  return pol;
}

Polynomial* PolynomialCopy(const Polynomial* pol)
{
  uint i;
  Polynomial *_cp = (Polynomial*) malloc(sizeof(Polynomial));
   if (!pol)
    return NULL;
   
   _cp->m_data = (double *) malloc((pol->m_degree+1) * sizeof(double));
   if (!_cp)
     return NULL; // tem de ter tratamento caso ocorra erro em uma situacao como essa, utilizar assert
  
  _cp->m_degree = pol->m_degree; 
  
  for(i=0; i < pol->m_degree+1; i++)
    _cp->m_data[i] = pol->m_data[i];
  
  return _cp;
}

Polynomial* PolynomialDelete(Polynomial* pol)
{
  if (!pol)
    return NULL;
  
  free(pol->m_data);
  pol = NULL;
  return;
}

double PolynomialFx(const Polynomial* pol, const double x)
{
  uint i;
  double results, data=0;
  
  for(i=0; i < pol->m_degree; i++){
    results =  pow(x,pol->m_degree-i);
    data+= (pol->m_data[i] *results);
  }
  
  data+= pol->m_data[i];
  return data;
}


void PolynomialShow(const Polynomial* pol)
{
  ulong i;
  
  printf("Polynomial degree %lu\n", pol->m_degree);
  
  for(i=0; i < pol->m_degree; i++){
    if (i)
      (pol->m_data[i] > 0) ? printf("+") : 0;
    
    printf(" %lf^%d ", pol->m_data[i], pol->m_degree-i);
  }

  (pol->m_data[i] > 0) ? printf("+") : 0;
  printf(" %lf\n", pol->m_data[i]);
  
  return;
}

void PolynomialSetConstants(Polynomial* pol, ...)
{
  va_list argPtr;
  va_start(argPtr, 0);
  uint i;
  
  if (!pol)
    return;
  
  for(i=0; i < (pol->m_degree+1); i++)
    pol->m_data[i] = va_arg(argPtr, double);
  
  va_end(argPtr);
  return;
}

void PolynomialSetConstantsAndDegree(Polynomial* pol, uint degree, ... )
{
  va_list argPtr;
  uint i;
  va_start(argPtr,0);
  if (!pol)
    return;
  
  pol->m_degree = degree;
  free(pol->m_data);
  
  for(i=0; i < (degree+1); i++)
   pol->m_data[i] = va_arg(argPtr, double);
  
  va_end(argPtr);
  return;
}