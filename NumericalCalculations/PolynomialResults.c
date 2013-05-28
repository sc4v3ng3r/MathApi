#include "PolynomialResults.h"

#define DEFAULT_SIZE 10

void PolynomialResultsAddData(PolynomialResultsTable* table, const double an, const double bn, 
			      const double xn, const double fxn, const double error, const uint iteration)
{
  
  if (!table)
    return;
  
  if (table->m_runner == (table->m_results+table->m_size) ){
    table->m_results = realloc(table->m_results, DEFAULT_SIZE * sizeof(ResultsInfo));
    table->m_size+=10;
  }
  
  table->m_total+=1;
  table->m_runner->m_interval.m_x = an;
  table->m_runner->m_interval.m_y = bn;
  table->m_runner->m_xn = xn;
  table->m_runner->m_fxn = fxn;
  table->m_runner->m_error = error;
  table->m_runner->m_iteration = iteration;
  table->m_runner++;
  return;
}

PolynomialResultsTable* PolynomialResultsTableInit()
{
  PolynomialResultsTable *table = (PolynomialResultsTable*) malloc(sizeof(PolynomialResultsTable));
  if (!table)
    return NULL;
  
  table->m_total = 0;
  table->m_results = (ResultsInfo *) malloc(DEFAULT_SIZE *sizeof(ResultsInfo));
  
  if (!table->m_results){
    free(table);
    return NULL;
  }
  
  table->m_runner = table->m_results;
  table->m_size = 10;
  return table;
}

void PolynomialResultsTableShow(const PolynomialResultsTable* table)
{
  ulong i;
  
  for(i=0; i < table->m_total; i++)
    printf("%lu\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t\n", table->m_results[i].m_iteration,
	   table->m_results[i].m_interval.m_x, table->m_results[i].m_interval.m_y,
	   table->m_results[i].m_xn, table->m_results[i].m_fxn, table->m_results[i].m_error);

  return;
}

void PolynomialResultsTableDelete(PolynomialResultsTable* table)
{
  free(table->m_results);
  table->m_results = NULL;
  table->m_runner = NULL;
  free(table);
  table = NULL;
}
