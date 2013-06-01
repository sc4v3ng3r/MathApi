#include "LinearSystem.h"

static ulong canChange(const unsigned long line, const unsigned long column, const Matrix *matrix);
static void trocarLinha(const ulong a, const ulong b, Matrix* matrix);
static int transformarEmZero(const ulong line, const ulong pivotRow, /*const int column,*/ double k, Matrix* matrix);
static int transformarEmUm(const ulong line, double k, Matrix* matrix );
static void square(Matrix* matrix);
static void rectangular(Matrix* matrix);
static uint checkDiagonal(Matrix *matrix);
static void findError(const ulong line, const ulong column, Matrix* matrix, double* previous);

LinearSystem* LinearSystemInit(const Matrix* matrix, BOOL vectorIndependent)
{
  ulong i;
  
  LinearSystem *linearSystem = (LinearSystem *) malloc(sizeof(LinearSystem));
  if (!linearSystem)
    return NULL;
  
  // verificador de erro!
  if (vectorIndependent){
    if (!matrix->m_type){ // if matrix is square with independent term vector
      linearSystem->m_systemMatrix = MatrixCopy(matrix);
      linearSystem->m_systemMatrix->m_type = RECTANGULAR;
      //linearSystem->m_systemMatrix->m_columns--;
    } 
    else {
      linearSystem->m_systemMatrix = MatrixCopy(matrix);
      linearSystem->m_systemMatrix->m_type = SQUARE;
    }
  }
  
  else {
    if (!matrix->m_type){ // if matrix is square and not have the independent term vector
      linearSystem->m_systemMatrix = MatrixCopy(matrix);
      
      for(i=0; i < matrix->m_rows; i++)
	linearSystem->m_systemMatrix->m_data[i] = realloc(linearSystem->m_systemMatrix->m_data[i],
							  ((matrix->m_columns+1) * sizeof(double)));
	linearSystem->m_systemMatrix->m_columns+=1;
    }
    else {
      linearSystem->m_systemMatrix = MatrixCopy(matrix);
      linearSystem->m_systemMatrix->m_type = SQUARE;
      
      for(i=0; i < matrix->m_rows; i++)
	linearSystem->m_systemMatrix->m_data[i] = realloc(linearSystem->m_systemMatrix->m_data[i],
							  ((matrix->m_columns+1) * sizeof(double)));
	linearSystem->m_systemMatrix->m_columns+=1;
    }
  }
  
  linearSystem->m_ca = linearSystem->m_cv = 0;
  linearSystem->m_solutionMatrix = NULL;
  return linearSystem;
}

LinearSystem* LinearSystemCopy(const LinearSystem* linearSystem)
{
  LinearSystem *_cp = (LinearSystem *) malloc(sizeof(LinearSystem));
  if (!_cp)
    return NULL;
  
  _cp->m_solutionMatrix = MatrixCopy(linearSystem->m_solutionMatrix);
  _cp->m_systemMatrix = MatrixCopy(linearSystem->m_systemMatrix);
  _cp->m_ca = linearSystem->m_ca;
  _cp->m_cv = linearSystem->m_cv;
  return _cp;
}

void LinearSystemDelete(LinearSystem* linearSystem)
{
  if (linearSystem->m_systemMatrix)
    MatrixDelete(linearSystem->m_systemMatrix);
  
  if (linearSystem->m_solutionMatrix)
    MatrixDelete(linearSystem->m_solutionMatrix);
  
  free(linearSystem);
  linearSystem = NULL;
  return;
}

void LinearSystemSetIndependentTermsVector(LinearSystem* linearSystem, ...)
{
  uint i;
  
  va_list argPtr;
  va_start(argPtr, 0);
  
  for(i=0; i <linearSystem->m_systemMatrix->m_rows; i++)
    linearSystem->m_systemMatrix->m_data[i][linearSystem->m_systemMatrix->m_columns-1] = va_arg(argPtr,double);
  
  va_end(argPtr);
  return;
}

void LinearSystemGaussJordan(LinearSystem* linearSystem)
{
  if (linearSystem->m_solutionMatrix)
    MatrixDelete(linearSystem->m_solutionMatrix);
  
  linearSystem->m_solutionMatrix = MatrixCopy(linearSystem->m_systemMatrix);
  switch(linearSystem->m_systemMatrix->m_type){
    case SQUARE:
      square(linearSystem->m_solutionMatrix);
      break;
    
    case RECTANGULAR:
      rectangular(linearSystem->m_solutionMatrix);
      break;
  }
  return;
}

InteractiveResultsTable* LinearSystemJacobi(const LinearSystem* linearSystem, const ulong interations,
					   const double precision, const double* kick)
{
  ulong i,j,k; // contadores de loops
  ulong columnCount, rowCount;
  InteractiveResultsTable *results;
  double *beginKick, data; // errMIn is presicion
  
  columnCount = rowCount = 0;
  
  if (linearSystem->m_systemMatrix->m_type!= SQUARE) // se for matriz retangular ... cai fora!
    return NULL;
 
  if (!checkDiagonal(linearSystem->m_systemMatrix) ){
    puts("\n\t\tImpossivel Calcular esse sistema!\n");
    return NULL;
  }
  
  results = InteractiveResultsTableInit(interations,(linearSystem->m_systemMatrix->m_columns),TRUE);
   if (!results){
    printf("Error in JACOBI() results malloc [DEBUG] %s\n", strerror(errno));
    return;
  }
  
  results->m_precsion = precision;
  
  beginKick = (double*) malloc(linearSystem->m_systemMatrix->m_rows*sizeof(double));
  if (!beginKick){
    printf("Error in JACOBI() beginKick malloc [DEBUG] %s\n", strerror(errno));
    return;
  }
  
  for(i=0; i < linearSystem->m_systemMatrix->m_rows; i++)
    beginKick[i] = kick[i];
  
   rowCount=0;
   for(k=0; k < interations; k++){ // Loop de iteracoes
    results->m_total+=1;
    columnCount=0;
    
    for(i=0; i<linearSystem->m_systemMatrix->m_rows;i++){ // percorre as linhas (equacoes)
      double sum = linearSystem->m_systemMatrix->m_data[i][linearSystem->m_systemMatrix->m_columns-1]; // obtenho o termo independente
      
//       printf("Termo independente Linha %lu : %lf\n", i, sum);
      for(j=0;j<linearSystem->m_systemMatrix->m_columns;j++){ // percorre as colunas [LOOP DE SOMATORIO!]
	if(j==i)
	  continue;
	sum+= ( linearSystem->m_systemMatrix->m_data[i][j]* -1 ) * beginKick[j];
      }
//       printf("Fazendo %lf / %lf\n", sum,linearSystem->m_systemMatrix->m_data[i][i]);
      results->m_results->m_data[rowCount][columnCount] = sum/linearSystem->m_systemMatrix->m_data[i][i];
      columnCount+=1;
    }
    
    findError(rowCount,columnCount-1,results->m_results, beginKick); // encontra o erro, e ja coloca na matriz!
    
    // atualizando os valores que serao utilizados na proxima iteracao!
    for(i=0; i < linearSystem->m_systemMatrix->m_rows;i++)
	beginKick[i] = results->m_results->m_data[rowCount][i];
    
    // verificar a convergencia! e caso seja < que a tolerancia, sair do loop
     // colocando os novos valores para serem utilizados na proxima iteracao
    data = results->m_results->m_data[rowCount][(linearSystem->m_systemMatrix->m_columns-1)];
    
    if (data < precision)
      break;
      
     rowCount+=1;
   }//final do loop de iteracoes  
   
  free(beginKick);
  results->m_operation = JACOBI;
  return results;
}

InteractiveResultsTable *LinearSystemGaussSeidel(const LinearSystem* linearSystem, const ulong interations,
						const double precision, const double* kick)
{
  ulong i,j,k; // contadores de loops
  ulong columnCount, rowCount;
  InteractiveResultsTable* results;
  double* beginKick, *kickBackup, data; // errMIn is presicion
  
  columnCount = rowCount = 0;// talvez nao necessario!
  
  if (!linearSystem)
    return NULL;
  
  if (linearSystem->m_systemMatrix->m_type!= SQUARE) // se for matriz retangular ... cai fora!
    return NULL;
 
   // Checar, pois nao funciona direito! WARNING
  if ( !checkDiagonal(linearSystem->m_systemMatrix) ){ // VAI MUDAR PARA O CRITERIO DE CONVERGENCIA
    puts("\n\t\tImpossivel Calcular esse sistema!\n");
    return NULL;
  }
  
  results = InteractiveResultsTableInit(interations,linearSystem->m_systemMatrix->m_columns,TRUE);
  if (!results){
    printf("InteractiveResultsTable LinearSystemGaussSeidel ERROR %s\n", strerror(errno));
    return NULL;
  }
  
  beginKick = (double*) malloc(linearSystem->m_systemMatrix->m_rows*sizeof(double));
  
  for(i=0; i < linearSystem->m_systemMatrix->m_rows; i++)
    beginKick[i] = kick[i];
  
  kickBackup = (double*) malloc(linearSystem->m_systemMatrix->m_rows*sizeof(double));
  
  for(i=0; i<linearSystem->m_systemMatrix->m_rows; i++)
    kickBackup[i] = beginKick[i];
  
  rowCount=0;
  results->m_operation = GAUSS_SEIDEL;
  results->m_precsion = precision;
  
  for(k=0; k < interations; k++){ // Loop de iteracoes
    columnCount=0;
    results->m_total+=1;
    for(i=0; i<linearSystem->m_systemMatrix->m_rows;i++){ // percorre as linhas (equacoes)
      double sum = linearSystem->m_systemMatrix->m_data[i][linearSystem->m_systemMatrix->m_columns-1]; // obtenho o termo independente
      for(j=0;j<linearSystem->m_systemMatrix->m_columns;j++){ // percorre as colunas [LOOP DE SOMATORIO!
	if(j==i)
	  continue;
	sum+= ( linearSystem->m_systemMatrix->m_data[i][j]*-1 ) *kickBackup[j];
      }
      kickBackup[i] = sum/linearSystem->m_systemMatrix->m_data[i][i];
      results->m_results->m_data[rowCount][columnCount] = kickBackup[i];
      columnCount+=1;
    }
    
    findError(rowCount,columnCount-1,results->m_results, beginKick); // encontra o erro, e ja coloca na matriz!
    
    for(i=0; i < linearSystem->m_systemMatrix->m_rows;i++)
      beginKick[i] = kickBackup[i];
    
    // verificar a convergencia! e caso seja < que a tolerancia, sair do loop
     // colocando os novos valores para serem utilizados na proxima iteracao
    data = results->m_results->m_data[rowCount][(results->m_results->m_columns-1)];
    if ( data < precision )
      break;
    rowCount+=1;
   }//final do loop de iteracoes
   
   free(beginKick);
   free(kickBackup);
   return results;
}

// trocar de linha vai ter de verificar se foi possivel trocar, e retornar verdadeiro se sim, ou false se nao!
static void trocarLinha(const ulong a, const ulong b, Matrix* matrix)
{
  unsigned long i;
  double *backup = (double*) malloc(matrix->m_columns* sizeof(double));
  
  for(i=0; i < matrix->m_columns; i++){
    backup[i] = matrix->m_data[a][i];
    matrix->m_data[a][i] = matrix->m_data[b][i];
  }
  
  for(i=0; i < matrix->m_columns; i++)
    matrix->m_data[b][i] = backup[i];
  
  free(backup);
  backup = NULL;
  return;
}
static int transformarEmUm(const ulong line, double k, Matrix* matrix)
{
  ulong j;
  if (!k)
    return FALSE;
  
  for(j=0; j < matrix->m_columns; j++){
    matrix->m_data[line][j]/=k;
     if (matrix->m_data[line][j] == (-0) )
      matrix->m_data[line][j]=0;
  }
  return TRUE;
}
// linha onde vai comecar a zerar, linha pivo, constante, e a matriz
static int transformarEmZero(const ulong line, const ulong pivotRow, /*const int column,*/ double k, Matrix* matrix)
{
  ulong i;
  double data;
  
  if (!k)
    return FALSE;
  
  k*=-1;

  for(i=0; i < matrix->m_columns; i++){
    data = k* matrix->m_data[pivotRow][i]; 
    matrix->m_data[line][i]+=data;
    if (matrix->m_data[line][i] == (-0) )
      matrix->m_data[line][i]=0;
  }
  return TRUE;
}

static void square(Matrix* matrix)
{
  ulong i, rowCount, columnCount;
  rowCount = columnCount = 0;
  
  while( (rowCount < matrix->m_rows) && (columnCount < (matrix->m_columns-1)) ){
    if (matrix->m_data[rowCount][columnCount]){
      transformarEmUm(rowCount,matrix->m_data[rowCount][columnCount],matrix);
      for(i=0; i < matrix->m_rows;i++){
	if (i==rowCount)
	  continue;
	transformarEmZero(i,rowCount/*pivo!*/,matrix->m_data[i][columnCount],matrix);
      }
    }
    // se tudo acontecer normal o primeiro passo acaba aki!
    else { // caso o valor seja zero, procurar alguem para trocar, caso nao encontre, segue adiante
      ulong line;
      for(i=0; i < matrix->m_rows; i++){
	line = canChange(i,columnCount,matrix);
	if (line){
	  trocarLinha(line,rowCount,matrix);
	  transformarEmUm(rowCount,matrix->m_data[rowCount][columnCount],matrix);
	  for(i=0; i < matrix->m_rows;i++){
	    if (i==rowCount)
	      continue;
	    transformarEmZero(i,rowCount/*pivo!*/,matrix->m_data[i][columnCount],matrix);
	  }
	  break;
	}
      }
    }
    rowCount++;
    columnCount++;
  }
}

void rectangular(Matrix* matrix)
{
  ulong i, rowCount, columnCount;
  rowCount = columnCount = 0;
 // puts("Matrix Retangular");
  while ( (rowCount < matrix->m_rows) && (columnCount < (matrix->m_columns-1) )){
    if (matrix->m_data[rowCount][columnCount]){
      transformarEmUm(rowCount,matrix->m_data[rowCount][columnCount],matrix);
      for(i=0; i < matrix->m_rows;i++){
	if (i==rowCount)
	  continue;
	transformarEmZero(i,rowCount/*pivo!*/,matrix->m_data[i][columnCount],matrix);
      }
      rowCount++;
      columnCount++;
    } // se tudo acontecer normal o primeiro passo acaba aki!
    else { // caso o valor seja zero, procurar alguem para trocar, caso nao encontre, segue adiante
      ulong line;
      for(i=0; i < matrix->m_rows; i++){
	line = canChange(i,columnCount,matrix);
	if (line){
	  trocarLinha(line,rowCount,matrix);
	  transformarEmUm(rowCount,matrix->m_data[rowCount][columnCount],matrix);
	  for(i=0; i < matrix->m_rows;i++){
	    if (i==rowCount)
	      continue;
	    transformarEmZero(i,rowCount/*pivo!*/,matrix->m_data[i][columnCount],matrix);
	  }
	  rowCount++;
	  columnCount++;
	  break;
	}
      }
      if (!line)
	columnCount++;
    }
  }
}

static ulong canChange(const unsigned long line, const unsigned long column, const Matrix *matrix){
  
  ulong i= line+1;
  if (i >= matrix->m_rows)
    return FALSE;
  
  while(i < matrix->m_rows){
    if (matrix->m_data[i][column])
      return i;
    i+=1;
  }
  return FALSE;
}

static uint checkDiagonal(Matrix *matrix){
    uint i;
    
    for(i=0; i< matrix->m_rows; i++){
      if (!matrix->m_data[i][i])
	trocarLinha(i,i+1,matrix);
    }
    
    for(i=0; i <matrix->m_rows;i++)
      if (!matrix->m_data[i][i])
	return FALSE;
      
  return TRUE;
}

static void findError(const ulong line, const ulong column, Matrix* matrix, double* previous){
  ulong i;
  double* erros, data;
  
  erros = (double*) malloc((matrix->m_columns-1) * sizeof(double));
  
  for(i=0; i < matrix->m_columns-1;i++){
    erros[i] = matrix->m_data[line][i] - previous[i];

    if (erros[i] < 0)
      erros[i]*=-1;
  }
  // encontrando o maior erro! pra armazenar na matriz!
  data = erros[0];
  for(i=1; i < matrix->m_columns-1; i++){
    if (data < erros[i])
      data = erros[i];
  }
  
  matrix->m_data[line][column+1] = data;
  free(erros);
  return;
}