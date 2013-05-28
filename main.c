#include "./LinearAlgebra/LinearSystem.h"
#include "./MathBase/Polynomial.h"
#include "./NumericalCalculations/PolynomialRoot.h"
#include "./NumericalCalculations/PolynomialResults.h"

int main(){
  ulong i,j;
  
  LinearSystem *system;
  Matrix * matrix = MatrixInit(3,4);
  
  /*
  for(i=0; i < matrix->m_rows; i++)
    for(j=0; j < matrix->m_columns; j++)
      matrix->m_data[i][j] = (3.14 * (i+(j+1)) * 2) / 10;
    */
  
  matrix->m_data[0][0] = 1;
  matrix->m_data[0][1] = 1;
  matrix->m_data[0][2] = 2;
  matrix->m_data[0][3] = 9;
  
  matrix->m_data[1][0] = 0;
  matrix->m_data[1][1] = 2;
  matrix->m_data[1][2] = -7;
  matrix->m_data[1][3] = -17;
  
  matrix->m_data[2][0] = 3;
  matrix->m_data[2][1] = 6;
  matrix->m_data[2][2] = -5;
  matrix->m_data[2][3] = 0;
  
  
  MatrixShow(matrix);
  system = LinearSystemInit(matrix, TRUE);
  
  MatrixShow(system->m_systemMatrix);
  MatrixDelete(matrix);
  
  //LinearSystemSetIndependentTermsVector(system, 0.1,0.2,0.3,0.5,0.4,0.03);
  MatrixShow(system->m_systemMatrix);
  LinearSystemGaussJordan(system);
  MatrixShow(system->m_solutionMatrix);
  
  LinearSystemDelete(system);
  
  Polynomial *pol = PolynomialInit(3);
  PolynomialSetConstants(pol, 1.0, 0,-9.0, 3.0);
  OrderedPair *pair = (OrderedPair*) malloc(sizeof(OrderedPair));
  pair->m_x = 0.0;
  pair->m_y = 1.0;
  
  printf("Intervals %lf %lf\n", pair->m_x, pair->m_y);
  PolynomialResultsTable *table = PolynomialRootBissection(pol,pair,10,0.002);
  
  if (table){
    printf("%lu %lu\n", table->m_size, table->m_total);//  PolynomialResultsTableShow(table);
    PolynomialResultsTableDelete(table);
  }
  
  PolynomialShow(pol);
  PolynomialDelete(pol);
  free(pair);
  return 0;
}
