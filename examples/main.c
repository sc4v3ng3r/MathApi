#include "../src/Subjects/LinearAlgebra/LinearSystem.h"
#include "../src/Core/MathBase/Polynomial.h"
#include "../src/Subjects/NumericalCalculations/PolynomialRoot.h"
#include "../src/Subjects/NumericalCalculations/PolynomialResults.h"

/* this method main is only for test
 * 
 */

int main(){
  ulong i,j;
  
  LinearSystem *system;
  Matrix * matrix = MatrixInit(3,4);
  
  /*
  for(i=0; i < matrix->m_rows; i++)
    for(j=0; j < matrix->m_columns; j++)
      matrix->m_data[i][j] = (3.14 * (i+(j+1)) * 2) / 10;
    */
  // ONLY TEST ---------
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
//   -----------------------
  
  MatrixShow(matrix);
  system = LinearSystemInit(matrix, TRUE);
  
  MatrixShow(system->m_systemMatrix);
  MatrixDelete(matrix);
  
  //LinearSystemSetIndependentTermsVector(system, 0.1,0.2,0.3,0.5,0.4,0.03);
  MatrixShow(system->m_systemMatrix);
  LinearSystemGaussJordan(system);
  MatrixShow(system->m_solutionMatrix);
  LinearSystemDelete(system);
  
  Polynomial *pol = PolynomialInit(2);
  PolynomialSetConstants(pol, 1.0, 0.0,-3.0);
  
  OrderedPair *pair = (OrderedPair*) malloc(sizeof(OrderedPair));
  pair->m_x = 1;
  pair->m_y = 2;
  
  printf("Intervals %lf %lf\n", pair->m_x, pair->m_y);
  PolynomialResultsTable *table = PolynomialRootBissection(pol,pair,10,0.01);
  
  //printf("%lu %lf\n",table->m_results->m_iterator, table->m_results[1].m_data[1]);
  if (table){
    puts("\n Bissecao \n\n");
    PolynomialResultsTableShow(table);
    PolynomialResultsTableDelete(table);
  }
  
  Polynomial *pol2 = PolynomialInit(3);
  PolynomialSetConstants(pol2,1.0,0.0,(-9.0),3.0);
  
  PolynomialShow(pol2);
  pair->m_x = 0;
  pair->m_y = 1;
  
  table = PolynomialRootSecant(pol2, pair , 10, 0.0005);
  puts("\n SECANTE \n\n");
  PolynomialResultsTableShow(table);
  printf("\nRoot founded is: %lf\n", *table->m_root);
  
  PolynomialResultsTableDelete(table);
  PolynomialShow(pol);
  PolynomialDelete(pol);
  free(pair);
  PolynomialDelete(pol2);
  return 0;
}