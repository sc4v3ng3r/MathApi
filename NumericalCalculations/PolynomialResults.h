#ifndef POLYNOMIAL_RESULTS_H
#define POLYNOMIAL_RESULTS_H

#include "../MathBase/defines.h"
#include "../MathBase/OrderedPair.h"

typedef struct ResultsInfo {
  ulong m_iteration;
  OrderedPair m_interval;
  double m_xn, m_fxn, m_error;
}ResultsInfo;

typedef struct PolynomialResultsTable {
  ulong m_total, m_size;
  ResultsInfo * m_results;
  ResultsInfo * m_runner;
}PolynomialResultsTable;

PolynomialResultsTable * PolynomialResultsTableInit();
void PolynomialResultsTableShow(const PolynomialResultsTable* table);
void PolynomialResultsTableDelete(PolynomialResultsTable* table);

void PolynomialResultsAddData(PolynomialResultsTable *table, const double an, const double bn,
			      const double xn, const double fxn, const double error, const uint iteration);
#endif