#include "InteractiveResultsTable.h"

InteractiveResultsTable* InteractiveResultsTableInit(const ulong rows, const ulong columns, BOOL error)
{
  ushort tester = 0;
  
  InteractiveResultsTable *table = (InteractiveResultsTable*) malloc(sizeof(InteractiveResultsTable));
  if (!table){
    printf("InteractiveResultsTableInit ERROR %s\n", strerror(errno));
    return NULL;
  }
  if (!error)
    tester+=1;
  
  table->m_results = MatrixInit(rows, (columns/*+tester*/));
  table->m_operation = NOT_CALCULATED;
  table->m_precsion = 0;
  table->m_total = 0;
  table->m_solution;
  return table;
}

InteractiveResultsTable* InteractiveResultsTableCopy(const InteractiveResultsTable* table)
{
 
 return NULL;
}

void InteractiveResultsTableShow(const InteractiveResultsTable* table)
{
  ulong i;
 // Maybe here this switch is not nedded
  switch(table->m_operation){
    case GAUSS_SEIDEL:
    case JACOBI:
      for(i=0; i < table->m_results->m_columns-1; i++)
	printf("     X%lu \t", i);
  
      printf("  E\n");
      MatrixShowLines(table->m_results, table->m_total);
      printf("\n Precision: %lf\n", table->m_precsion);
      break;
      
    default:
      break;
  }
  return;
}

void InteractiveResultsTableDelete(InteractiveResultsTable* table)
{
  if(!table)
    return;
  
  MatrixDelete(table->m_results);
  free(table);
  return;
}
