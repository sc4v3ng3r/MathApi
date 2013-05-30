#ifndef RESULTS_TABLE_H
#define RESULTS_TABLE_H

#include "defines.h"
#include "OrderedPair.h"

typedef enum Operation {JACOBI=1, BISSECTION=2, SECANT=3 }Operation;

typedef struct ResultsTable{
  uint m_iterator;
  OrderedPair m_pair;
  double *m_data;
  double m_error;
  Operation m_lastOperation;
}ResultsTable;

ResultsTable * ResultsTableInit();
void ResultsTableDelete(ResultsTable* table);
void ResultsTableShow(const ResultsTable* table);
void ResultsTableAddData(ResultsTable* table, const uint iterator, const OrderedPair pair,
			 const double* data, const double error, const Operation operation);
#endif