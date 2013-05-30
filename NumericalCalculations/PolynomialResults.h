#ifndef POLYNOMIAL_RESULTS_H
#define POLYNOMIAL_RESULTS_H

#include "../MathBase/defines.h"
#include "../MathBase/OrderedPair.h"
#include "../MathBase/ResultsTable.h"
/*
typedef struct ResultsInfo {
  ulong m_iteration;
  OrderedPair m_interval;
  double m_xn, m_fxn, m_error; // maybe here has to be a double pointer, and we need of one metaStruct...
}ResultsInfo;
*/

typedef struct PolynomialResultsTable {
  ulong m_total, m_size;
  ResultsTable * m_results;
  ResultsTable * m_runner;
  double *m_root; // this is a shortcut for the last root founded!
}PolynomialResultsTable;

PolynomialResultsTable * PolynomialResultsTableInit(const ulong maxResults);
void PolynomialResultsTableShow(const PolynomialResultsTable* table);
void PolynomialResultsTableDelete(PolynomialResultsTable* table);
void PolynomialResultsAddData(PolynomialResultsTable *table, const OrderedPair pair,const double* data,
			      const double error, const uint iteration, Operation op);
#endif