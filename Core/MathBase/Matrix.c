#include "Matrix.h"

Matrix* MatrixInit(const ulong rows, const ulong columns)
{ 
  ulong i;
  Matrix *data=0;
  if(! (data = (Matrix*) malloc(1*sizeof(Matrix))) ){
    printf("Memory alocation [Matrix*] failed in | DEBUG[INIT] %s= \n", strerror(errno));
    return NULL;
  }
  
  (rows == columns) ? (data->m_type = SQUARE) : (data->m_type = RECTANGULAR);
  
  data->m_rows = rows;
  data->m_columns= columns;
  
  if(! (data->m_data = (double**) malloc(rows*sizeof(double*))) ){
    printf("Memory alocation [ROWS]failed in | DEBUG[INIT] double %s\n", strerror(errno));
    return NULL;
  }
  
  for(i=0; i < rows; i++)
    data->m_data[i] = (double*) malloc(columns*sizeof(double));    
  
  MatrixReset(data);
  return data;
}

void MatrixDelete(Matrix* matrix)
{
  ulong i;
  for(i=0; i < matrix->m_rows; i++)
    free(matrix->m_data[i]);
  
  free(matrix->m_data);
  free(matrix);
  matrix = NULL;
  return;
}
Matrix* MatrixCopy(const Matrix* matrix)
{
  ulong i,j;
  Matrix* _cp; // copy
  
  if (matrix == NULL)
    return NULL;
  
  _cp = MatrixInit(matrix->m_rows,matrix->m_columns);
  if (!_cp)
    return NULL;
  
  for(i=0;i < matrix->m_rows; i++)
    for(j=0; j < matrix->m_columns; j++)
      _cp->m_data[i][j] = matrix->m_data[i][j];
  
  return _cp;
}

void MatrixReset(Matrix* matrix)
{
  ulong i,j;
  
  if (matrix)
    for(i=0; i < matrix->m_rows; i++)
      for(j=0; j < matrix->m_columns; j++)
	matrix->m_data[i][j] = 0x0;
  
  return;
}

void MatrixShow(const Matrix* matrix)
{
  ulong i,j;
  
  for(i=0; i< matrix->m_rows; i++){
    printf("L %lu ",i);
    for(j=0; j< matrix->m_columns; j++)
      printf("%.10lf  ", matrix->m_data[i][j]);
    puts("");
  }
  puts("\n");
  return;
}
