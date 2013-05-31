#ifndef POLYNOMIAL_RESULTS_H
#define POLYNOMIAL_RESULTS_H

#include "../../Core/MathBase/defines.h"
#include "../../Core/MathBase/OrderedPair.h"
#include "../../Core/ResultsTables/ResultsTable.h"

typedef struct PolynomialResultsTable {
  ulong m_total, m_size;
  ResultsTable * m_results;
  ResultsTable * m_runner;
  double *m_root; // this is a shortcut for the last root founded!
}PolynomialResultsTable;

typedef PolynomialResultsTable prt;

PolynomialResultsTable * PolynomialResultsTableInit(const ulong maxResults);
void PolynomialResultsTableShow(const PolynomialResultsTable* table);
void PolynomialResultsTableDelete(PolynomialResultsTable* table);
void PolynomialResultsAddData(PolynomialResultsTable *table, const OrderedPair pair,const double* data,
			      const double error, const uint iteration, Operation op);
#endif