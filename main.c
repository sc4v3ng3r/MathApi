//#include "../../src/MathBase/Matrix.h"
#include "src/Core/ResultsTables/InteractiveResultsTable.h"
#include "src/Subjects/LinearAlgebra/LinearSystem.h"

int main(){
  
  Matrix *matrix = MatrixInit(2,3);
  double *data = (double*) malloc(2* sizeof(double));
  data[0] = 0;
  data[1] = 0;
  
  matrix->m_data[0][0] = 2;
  matrix->m_data[0][1] = -1;
  matrix->m_data[0][2] = 1;
  matrix->m_data[1][0] = 1;
  matrix->m_data[1][1] = 2;
  matrix->m_data[1][2] = 3;
  
  
  LinearSystem *linearSystem = LinearSystemInit(matrix,1);
  MatrixShow(linearSystem->m_systemMatrix);
  
  InteractiveResultsTable* table = LinearSystemJacobi(linearSystem,10,0.01, data);
  
  if (table){
    InteractiveResultsTableShow(table);
    InteractiveResultsTableDelete(table);
  } else puts("FUUUUUUCKKKK");
  LinearSystemDelete(linearSystem);
  MatrixDelete(matrix);
  return 0;
}