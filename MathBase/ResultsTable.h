#ifndef RESULTS_TABLE_H
#define RESULTS_TABLE_H

#include "defines.h"
#include "OrderedPair.h"

#define BISSECTION_PRINT printf("K\t    AN\t          BN\t           XN\t           F(XN)\t        E\n");

typedef enum Operation {JACOBI=1, BISSECTION=2, SECANT=3 }Operation;

typedef struct ResultsTable{
  uint m_iterator;
  OrderedPair m_pair;
  double *m_data;
  double m_error;
  double m_precision;
  Operation m_lastOperation;
}ResultsTable;

ResultsTable * ResultsTableInit();
void ResultsTableDelete(ResultsTable* table);
void ResultsTableShow(const ResultsTable* table);
void ResultsTableAddData(ResultsTable* table, const uint iterator, const OrderedPair pair,
			 const double* data, const double error, const Operation operation);
#endif