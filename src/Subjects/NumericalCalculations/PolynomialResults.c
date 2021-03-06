#include "PolynomialResults.h"

#define DEFAULT_SIZE 10

void PolynomialResultsAddData(PolynomialResultsTable* table,const OrderedPair pair,const double* data,
			      const double error, const uint iteration, Operation op)
{
  ulong i;
  
  if (!table)
    return;
  
  //printf("Table 0x%x Runner 0x%x  bound is : 0x%x\n", table->m_results, table->m_runner, (table->m_results+table->m_size));
  
  if (table->m_runner == (table->m_results+table->m_size) ){
    return;
    /*
    table->m_results = (ResultsTable*)realloc(table->m_results, 10*sizeof(ResultsTable));
    
    table->m_size+=10;
    for(i=table->m_total; i < table->m_size; i++)
      table->m_results[i].m_data = NULL;
    
    puts("New 10 alocation\n");*/
  } 
  ResultsTableAddData(table->m_runner,iteration,pair, data, error, op);
  
  table->m_total+=1;
  table->m_runner++;
  return;
}

PolynomialResultsTable* PolynomialResultsTableInit(const ulong maxResults)
{
  PolynomialResultsTable *table = (PolynomialResultsTable*) malloc(sizeof(PolynomialResultsTable));
  if (!table){
    printf("PolynomialResultsTableInit ERROR %s\n", strerror(errno));
    return NULL;
  }
  
  table->m_total = 0;
  //   WARNING ALTERNATIVE SOLUTION, REPAIR THE REALLOC PROBLEM IN ResultsTable
  //table->m_results = ResultsTableInit();
  
  table->m_results = (ResultsTable*) malloc(maxResults* sizeof(ResultsTable));
   if (!table->m_results){
    free(table);
    return NULL;
  }
  
  table->m_runner = table->m_results;
  ulong i; for(i=0; i < maxResults; i++) table->m_results[i].m_data = NULL;
  // END OF ALTERNATIVE SOLUTION
  
  table->m_size = maxResults;
  return table;
}

void PolynomialResultsTableShow(const PolynomialResultsTable* table)
{
  ulong i;
  
  switch(table->m_results->m_lastOperation){
    case BISSECTION:
      BISSECTION_PRINT
      for(i=0; i < table->m_total; i++)
	ResultsTableShow(&table->m_results[i]);
      printf("\nTotal iterations: %lu\nMax iterations: %lu\nRoot founded: %lf\n",
	     table->m_total, table->m_size,*table->m_root);
      break;
      
    case SECANT:
      SECANT_PRINT
      for(i=0; i < table->m_total; i++)
	ResultsTableShow(&table->m_results[i]);
      break;
  }
  return;
}

void PolynomialResultsTableDelete(PolynomialResultsTable* table)
{
  ulong i;
  
  for(i=0; i < table->m_size; i++){
   // printf("PolynomialResultsTableDelete clear 0x%x\n", &table->m_results[i]);
    free(table->m_results[i].m_data);
    table->m_results[i].m_data = NULL;
  }
  
  ResultsTableDelete(table->m_results);
  free(table);
  table = NULL;
}