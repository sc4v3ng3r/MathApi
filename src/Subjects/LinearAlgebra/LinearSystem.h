#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "../../Core/MathBase/defines.h"
#include "../../Core/MathBase/Matrix.h"
#include "../../Core/ResultsTables/ResultsTable.h"

#include <stdarg.h>

typedef struct LinearSystem {
  Matrix *m_systemMatrix;
  Matrix *m_solutionMatrix;
  ulong m_ca, m_cv;
}LinearSystem;

typedef struct InterativeResultsTable {
  ulong m_size;
  ulong m_total;
  OrderedPair m_pair;
  ResultsTable *m_retsults;
}InterativeResultsTable;

LinearSystem * LinearSystemInit(const Matrix* matrix, BOOL vectorIndependent);
LinearSystem * LinearSystemCopy(const LinearSystem* linearSystem);
void LinearSystemDelete(LinearSystem* linearSystem);
// iterators methods
void LinearSystemJacobi(LinearSystem* linearSystem);
void LinearSystemGaussSeidel(LinearSystem* linearSystem);
void LinearSystemSetIndependentTermsVector(LinearSystem* linearSystem, ...);

void LinearSystemGaussJordan(LinearSystem* linearSystem);
// LinearSystem* LinearSystemCreate(const ulong rows, const ulong columns);
#endif