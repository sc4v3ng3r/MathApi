#include "ResultsTable.h"
#define RESULTS_TABLE_ADD_ERROR printf("ResultsTableAddData ERROR %s\n", strerror(errno));

ResultsTable* ResultsTableInit()
{
  ulong i;
  ResultsTable *table = (ResultsTable*) malloc(sizeof(ResultsTable));
  if (!table){
    printf("RESULTS TABLE INIT ERROR %s\n", strerror(errno));
    return NULL;
  }
  table->m_data = NULL;
  table->m_dataSize=0;
  return table;
}

void ResultsTableShow(const ResultsTable* table)
{
  //uint i,j;
  if (!table){
    puts("Empty table!");
    return;
  }
  
  switch(table->m_lastOperation){
    case BISSECTION:
      printf("%lu\t %lf\t %lf\t %lf\t %lf\t %lf\n",table->m_iterator,table->m_pair.m_x, table->m_pair.m_y,
	     table->m_data[0], table->m_data[1], table->m_error); 
      break;
      
    case SECANT:
      printf("%lu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",table->m_iterator,table->m_pair.m_x,
	     table->m_pair.m_y, table->m_data[0], table->m_data[1],table->m_data[2],table->m_data[3]);
      break;
      
    default:
      break;
  }

  return;
}

void ResultsTableAddData(ResultsTable* table, const uint iterator,const OrderedPair pair,
			 const double* data, const double error, const Operation operation)
{
  table->m_pair.m_x = pair.m_x;
  table->m_pair.m_y = pair.m_y;
  table->m_iterator = iterator;
  table->m_error = error;
  table->m_lastOperation = operation;
  
  switch(operation){
    case BISSECTION:
      table->m_data = (double*) malloc(2*sizeof(double));
      if (!table->m_data){
	RESULTS_TABLE_ADD_ERROR
	break;
      }
      table->m_dataSize = 2;
      memcpy(table->m_data,data,2*sizeof(double));
      break;
      
    case JACOBI:
      break;
      
    case SECANT:
      table->m_data = (double*) malloc(4*sizeof(double));
      if (!table->m_data){
	RESULTS_TABLE_ADD_ERROR
	break;
      }
      table->m_dataSize = 4;
      memcpy(table->m_data, data, 4*sizeof(double));
      break;
  }
  return;
}

ResultsTable* ResultsTableCopy(const ResultsTable* table)
{
  ulong i;
  
  ResultsTable *cp = ResultsTableInit();
  if (!cp)
    return NULL;
  
  if (table->m_dataSize){
    cp->m_data = (double*) malloc(table->m_dataSize*sizeof(double));
    if (!cp->m_data)
      return NULL;
  }
  
  cp->m_iterator = table->m_iterator;
  cp->m_lastOperation = table->m_lastOperation;
  cp->m_pair = table->m_pair;
  cp->m_precision = table->m_precision;
  cp->m_error = table->m_error;
  cp->m_dataSize = table->m_dataSize;
  
  memcpy(cp->m_data, table->m_data, table->m_dataSize*sizeof(double));
  return cp;
}

void ResultsTableDelete(ResultsTable* table)
{
  
  if (table->m_data){
//     printf("Yes 0x%x data 0x%x is NOT NULL\n", table, table->m_data);
    free(table->m_data);
  }
  
  free(table);
  table = NULL;
  return;
}
