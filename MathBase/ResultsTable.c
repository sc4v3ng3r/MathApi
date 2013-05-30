#include "ResultsTable.h"

ResultsTable* ResultsTableInit()
{
  ulong i;
  ResultsTable *table = (ResultsTable*) malloc(sizeof(ResultsTable));
  if (!table){
    printf("RESULTS TABLE INIT ERROR %s\n", strerror(errno));
    return NULL;
  }
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
    
    case JACOBI:
      
      break;
      
    case BISSECTION:
      printf("%lu\t %lf\t %lf\t %lf\t %lf\t %lf\n",table->m_iterator,table->m_pair.m_x, table->m_pair.m_y,
	     table->m_data[0], table->m_data[1], table->m_error); 
      break;
      
    case SECANT:
      
      break;
      
    default:
      break;
  }

  return;
}

void ResultsTableAddData(ResultsTable* table, const uint iterator,const OrderedPair pair,
			 const double* data, const double error, const Operation operation)
{
  table->m_pair = pair;
  table->m_iterator = iterator;
  table->m_error = error;
  table->m_lastOperation = operation;
  
  switch(operation){
    case BISSECTION:
      table->m_data = (double*) malloc(2*sizeof(double));
      if (!table->m_data){
	printf("ResultsTableAddData ERROR %s\n", strerror(errno));
	break;
      }
      memcpy(table->m_data,data,2*sizeof(double));
      break;
      
    case JACOBI:
      break;
      
    case SECANT:
      break;
  }
  return;
}

void ResultsTableDelete(ResultsTable* table)
{
  
  if (table->m_data){
    printf("Yes 0x%x data 0x%x is NOT NULL\n", table, table->m_data);
    free(table->m_data);
  }
  
  free(table);
  table = NULL;
  return;
}
