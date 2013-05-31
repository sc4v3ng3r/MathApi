#ifndef POLYNOMIAL_ROOT_H
#define POLYNOMIAL_ROOT_H

#include "../../Core/MathBase/Polynomial.h"
#include "../../Core/MathBase/OrderedPair.h"
#include "PolynomialResults.h"

PolynomialResultsTable * PolynomialRootBissection(const Polynomial* pol, const OrderedPair *interval,
			      const uint iterators, const double precision);
PolynomialResultsTable * PolynomialRootSecant(const Polynomial* pol, const OrderedPair *interval,
					      const uint iterators, const double precision);

#endif