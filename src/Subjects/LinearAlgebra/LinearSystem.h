#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "../../Core/MathBase/defines.h"
#include "../../Core/MathBase/Matrix.h"
#include "../../Core/ResultsTables/ResultsTable.h"
#include "../../Core/ResultsTables/InteractiveResultsTable.h"

#include <stdarg.h>

typedef struct LinearSystem {
  Matrix *m_systemMatrix;
  Matrix *m_solutionMatrix;
  ulong m_ca, m_cv;
}LinearSystem;

LinearSystem * LinearSystemInit(const Matrix* matrix, BOOL vectorIndependent);
LinearSystem * LinearSystemCopy(const LinearSystem* linearSystem);
void LinearSystemDelete(LinearSystem* linearSystem);
// iteractive methods
InteractiveResultsTable* LinearSystemJacobi(const LinearSystem* linearSystem, const ulong interations,
					   const double precision, const double* kick);
InteractiveResultsTable *LinearSystemGaussSeidel(const LinearSystem* linearSystem, const ulong interations,
						const double precision, const double* kick);

void LinearSystemSetIndependentTermsVector(LinearSystem* linearSystem, ...);

void LinearSystemGaussJordan(LinearSystem* linearSystem);
// LinearSystem* LinearSystemCreate(const ulong rows, const ulong columns);
#endif