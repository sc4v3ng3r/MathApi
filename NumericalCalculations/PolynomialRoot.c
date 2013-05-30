#include "PolynomialRoot.h"

static BOOL PolynomialRootIntervalValidador(const Polynomial *pol, const OrderedPair* interval);

PolynomialResultsTable * PolynomialRootBissection(const Polynomial* pol, const OrderedPair* interval,
						  const uint iterators, const double precision)
{
  ulong i;
  double error, *data;
  OrderedPair pair;
  
  PolynomialResultsTable *table = PolynomialResultsTableInit(iterators);
  
  if (!table)
    return NULL;
  
  if (!PolynomialRootIntervalValidador(pol, interval)){
    printf("PolynomialRootBissection [ERROR] BAD INTERVAL %lf %lf\n", interval->m_x, interval->m_y);
    return NULL;
  }
  
  double err = (interval->m_y - interval->m_x);
  if (err < 0)
    err*=-1;
  
  if ( err < precision){
    printf("(b-a) < E\n");
    // ENcontrar qualquer X para EM [a,b];
  }
  
  pair.m_x = interval->m_x;
  pair.m_y = interval->m_y;
  
  data = (double *) malloc (2*sizeof(double));
  
  for(i=0; i < iterators; i++){
    
    data[0] = (pair.m_x + pair.m_y) / 2; // Xi
    data[1] = PolynomialFx(pol,data[0]); // F(xi)
    
    error = (data[0] - pair.m_y);
    if (error < 0)
      error*=-1;
       
    PolynomialResultsAddData(table, pair, data, error, i, BISSECTION);
    
    if (PolynomialFx(pol,pair.m_x) * data[1] > 0)
      pair.m_x = data[0];
    
    else pair.m_y = data[0];
    
    // final da sessao de testes de definicoes
    if (error < precision)
      break;
  }
  
  table->m_root = &table->m_results[table->m_total-1].m_data[0];
  return table;
  //SYSTEM_PAUSE
}

static BOOL PolynomialRootIntervalValidador(const Polynomial* pol, const OrderedPair* interval)
{
  
  if ( ( (PolynomialFx(pol,interval->m_x)) * (PolynomialFx(pol,interval->m_y)) ) < 0 )
    return TRUE;
  
  return FALSE;
}

// WARNING I AM WORKING HERE!
PolynomialResultsTable* PolynomialRootSecant(const Polynomial* pol, const OrderedPair* interval,
					     const uint iterators, const double precision)
{
  ulong i;
  double Xn, pt1, pt2, pt3, *data;
  OrderedPair pair;
  
  PolynomialResultsTable *table = PolynomialResultsTableInit(iterators);
  if (!table)
    return NULL;
  
  pair.m_x = interval->m_x;
  pair.m_y = interval->m_y;
  
  data = (double*) malloc(4*sizeof(double));
  if (!data)
    return NULL;
  
  for(i=0; i < iterators; i++){
    data[0] = PolynomialFx(pol,pair.m_x); // F(x0)
    data[1] = PolynomialFx(pol,pair.m_y); // F(x1)
    
    pt1 = pair.m_x * PolynomialFx(pol, pair.m_y);
    pt2 = pair.m_y * PolynomialFx(pol, pair.m_x);
    
    pt1 = pt1-pt2; // NUMERADOR
    
    pt3 = PolynomialFx(pol,pair.m_y) - PolynomialFx(pol, pair.m_x); // DENOMINADOR!
    Xn = pt1 / pt3;
    data[2] = Xn;
    data[3] = PolynomialFx(pol,Xn);
    
    PolynomialResultsAddData(table, pair, data, 0,i,SECANT);
    pt3 = PolynomialFx(pol,Xn);
    
    if (pt3 < 0)
      pt3*=-1;
    
    if ( pt3 < precision)
      break;
    
    pair.m_x = pair.m_y;
    pair.m_y = Xn;
  }
  table->m_root = &table->m_results[table->m_total-1].m_data[2];
  return table;
}