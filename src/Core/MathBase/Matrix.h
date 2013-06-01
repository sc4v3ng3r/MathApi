#ifndef MATRIX_H
#define MATRIX_H

#include "defines.h"

typedef enum Type { SQUARE=0, RECTANGULAR =1 } type;

typedef struct Matrix {
  double **m_data;
  ulong m_rows;
  ulong m_columns;
  type m_type;
}Matrix;

Matrix * MatrixInit(const ulong rows, const ulong columns);
void MatrixDelete(Matrix* matrix);
void MatrixShow(const Matrix* matrix);
void MatrixShowLines(const Matrix* matrix, const ulong totalLines);
Matrix * MatrixCopy(const Matrix* matrix);
void MatrixReset(Matrix* matrix);

#endif