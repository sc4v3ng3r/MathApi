#include "PolynomialRoot.h"

static BOOL PolynomialRootIntervalValidador(const Polynomial *pol, const OrderedPair* interval);

PolynomialResultsTable * PolynomialRootBissection(const Polynomial* pol, const OrderedPair* interval,
			      const uint iterators, const double precision)
{
  ulong i;
  double Xi, value, interA, interB, choosed, fa,fb,fxi;
  
  PolynomialResultsTable *table = PolynomialResultsTableInit();
  /*
  if (PolynomialRootIntervalValidador(pol, interval)){
    printf("PolynomialRootBissection [ERROR] BAD INTERVAL %lf %lf\n", interval->m_x, interval->m_y);
    return NULL;
  }*/
  
  printf("PolynomialRootBissection OK!\n");
  double err = (interval->m_y - interval->m_x);
  if (err < 0)
    err*=-1;
  
  if ( err < precision){
    printf("(b-a) < E\n");
    // ENcontrar qualquer X para EM [a,b];
  }
  
  //printf(" K\t    AN\t          BN\t           XN\t           F(XN)\t        E\n");
  interA = interval->m_x;
  interB = interval->m_y;
  
  for(i=0; i < iterators; i++){
    Xi = (interA + interB)/2;
    
    fa = PolynomialFx(pol,interA);
    fb = PolynomialFx(pol,interB);
    fxi = PolynomialFx(pol,Xi);
    
    // sessao de testes e definicoes
    value = (Xi - interB); // value is error!
    if (value < 0)
      value*=-1;
    
    printf("add | %lu %lf %lf %lf %lf %lf\n", i, interA,interB, Xi, fxi, value);
    PolynomialResultsAddData(table,interA,interB,Xi,fxi,value,i);
    
    if (fa * fxi > 0)
      interA = Xi;
    
    else interB = Xi;
    
    // final da sessao de testes de definicoes
    // iteracao,IterA, InterB, InterN, FX(IterN), (b-a)
    //printf(" %d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",i,interA,interB,Xi,fxi, value);
    if (value < precision)
      break;
  }
  
  return table;
  //SYSTEM_PAUSE
}

static BOOL PolynomialRootIntervalValidador(const Polynomial* pol, const OrderedPair* interval)
{
  if ( ( (PolynomialFx(pol,interval->m_x)) * (PolynomialFx(pol,interval->m_y)) ) < 0 )
    return 1;
  
  return 0;
}

PolynomialResultsTable* PolynomialRootSecant(const Polynomial* pol, const OrderedPair* interval,
					     const uint iterators, const double precision)
{
  PolynomialResultsTable *table;
  
  return table;
}
