#ifndef INTERACTIVE_RESULTS_TABLE_H
#define INTERACTIVE_RESULTS_TABLE_H

#include "ResultsTable.h"
#include "../MathBase/Matrix.h"

typedef struct InteractiveResultsTable {
  Matrix *m_results;
  ulong m_total;
  double m_precsion;
  double *m_solution;
  Operation m_operation;
}InteractiveResultsTable;

typedef InteractiveResultsTable irt;

InteractiveResultsTable *InteractiveResultsTableInit(const ulong rows, const ulong columns, BOOL error);
void InteractiveResultsTableShow(const InteractiveResultsTable* table);
InteractiveResultsTable * InteractiveResultsTableCopy(const InteractiveResultsTable* table);
void InteractiveResultsTableDelete(InteractiveResultsTable* table);
#endif