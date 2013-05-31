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

void LinearSystemJacobi(LinearSystem* linearSystem)
{
  int i,j,k; // loops counter
  ulong iterator, columnCount, rowCount;
  Matrix* results;
  double* beginKick, errMin, data; // errMIn is presicion
  columnCount = rowCount = 0;
  
  if (linearSystem->m_systemMatrix->m_type) // if matrix is not square, bail out!
    return;
 
   // Checar, pois nao funciona direito! WARNING
  if ( !checkDiagonal(linearSystem->m_systemMatrix) ){
    puts("\n\t\tImpossivel Calcular esse sistema!\n");
    return;
  }
  
  printf("Quantidade de iteracoes: ");
  scanf("%d",&iterator);
  
  results = MatrixInit(iterator,linearSystem->m_systemMatrix->m_columns-1);
   if (!results){
    printf("Error in JACOBI() results malloc [DEBUG] %s\n", strerror(errno));
    return;
  }
  
  beginKick = (double*) malloc(linearSystem->m_systemMatrix->m_rows*sizeof(double));
  if (!beginKick){
    printf("Error in JACOBI() beginKick malloc [DEBUG] %s\n", strerror(errno));
    return;
  }
  
  printf("\nTolerancia: ");
  scanf("%lf", &errMin);
  
  CLEAR_SCREEN
  
  for(i=0; i<linearSystem->m_systemMatrix->m_rows; i++){
    printf("Valor Inicial X[%d]: ",i);
    scanf("%lf", &beginKick[i]);
    CLEAR_SCREEN
  }
  
   rowCount=0;
   
   for(k=0; k < iterator; k++){ // Loop de iteracoes
    columnCount=0;
    
    for(i=0; i<linearSystem->m_systemMatrix->m_rows;i++){ // percorre as linhas (equacoes)
      double sum = linearSystem->m_systemMatrix->m_data[i][linearSystem->m_systemMatrix->m_columns-1]; // obtenho o termo independente
      
      for(j=0;j<linearSystem->m_systemMatrix->m_columns;j++){ // percorre as colunas [LOOP DE SOMATORIO!]
	if(j==i)
	  continue;
	sum+= ( linearSystem->m_systemMatrix->m_data[i][j]*-1 ) * beginKick[j];
      }
      
	results->m_data[rowCount][columnCount] = sum/linearSystem->m_systemMatrix->m_data[i][i];
	columnCount+=1;
    }
    
    findError(rowCount,columnCount-1,results, beginKick); // encontra o erro, e ja coloca na matriz!
    
    // atualizando os valores que serao utilizados na proxima iteracao!
    for(i=0; i < linearSystem->m_systemMatrix->m_rows;i++)
	beginKick[i] = results->m_data[rowCount][i];
    
    // verificar a convergencia! e caso seja < que a tolerancia, sair do loop
     // colocando os novos valores para serem utilizados na proxima iteracao
    
     data = results->m_data[rowCount][(results->m_columns-1)];
     
     printf("%.10f < %.10f \n", data,errMin);
     
     if ( data < errMin ){
       printf("Solucao encontrada!\n");
       break;
     }
      rowCount+=1;
   }//final do loop de iteracoes  
   
    CLEAR_SCREEN
    printf("Resultados das iteracoes:\n");
    
    //showMatrix(results);
    rowCount == iterator ? rowCount-=1, k-=1 : rowCount;
    printf("Solucao: ");
    
    for(i=0; i < results->m_columns-1; i++){
      printf("X[%d] = %.10f ", i, results->m_data[rowCount][i]);
    } 
    printf("Erro: %.10f\n Total de iteracoes: %d|\n", results->m_data[rowCount][i], k+1);
    (results->m_data[rowCount][i] < errMin) ? puts("Convergiu!") : puts("Nao convergiu!");
  
  MatrixDelete(results);
  free(beginKick);
  printf("Tolerancia %.10f\n", errMin);
  return;
}

void LinearSystemGaussSeidel(LinearSystem* linearSystem)
{
  ulong i,j,k; // contadores de loops
  ulong iterator, columnCount, rowCount;
  Matrix* results;
  double* beginKick, *kickBackup, errMin, data; // errMIn is presicion
  
  columnCount = rowCount = 0;// talvez nao necessario!
  
  if (linearSystem->m_systemMatrix->m_type) // se for matriz retangular ... cai fora!
    return;
 
   // Checar, pois nao funciona direito! WARNING
  if ( !checkDiagonal(linearSystem->m_systemMatrix) ){ // VAI MUDAR PARA O CRITERIO DE CONVERGENCIA
    puts("\n\t\tImpossivel Calcular esse sistema!\n");
    return;
  }
  
  printf("Quantidade de iteracoes: ");
  scanf("%d",&iterator);
  
  results = MatrixInit(iterator,linearSystem->m_systemMatrix->m_columns-1);
  
  beginKick = (double*) malloc(linearSystem->m_systemMatrix->m_rows*sizeof(double));
  kickBackup = (double*) malloc(linearSystem->m_systemMatrix->m_rows*sizeof(double));
  printf("\nTolerancia: ");
  scanf("%lf", &errMin);
  
  CLEAR_SCREEN
  
  // entrada para o chute inicial
  for(i=0; i<linearSystem->m_systemMatrix->m_rows; i++){
    printf("Valor Inicial X[%d]: ",i);
    scanf("%lf", &beginKick[i]);
    kickBackup[i] = beginKick[i];
    CLEAR_SCREEN
  }
  
   rowCount=0;
   
   for(k=0; k < iterator; k++){ // Loop de iteracoes
    columnCount=0;
    
    for(i=0; i<linearSystem->m_systemMatrix->m_rows;i++){ // percorre as linhas (equacoes)
      double sum = linearSystem->m_systemMatrix->m_data[i][linearSystem->m_systemMatrix->m_columns-1]; // obtenho o termo independente
      
      for(j=0;j<linearSystem->m_systemMatrix->m_columns;j++){ // percorre as colunas [LOOP DE SOMATORIO!]
	if(j==i)
	  continue;
	sum+= ( linearSystem->m_systemMatrix->m_data[i][j]*-1 ) *kickBackup[j];
      }
	kickBackup[i] = sum/linearSystem->m_systemMatrix->m_data[i][i];
	results->m_data[rowCount][columnCount] = kickBackup[i];
	// colocando os valores encontrados para a proxima iteracao... :)
	columnCount+=1;
    }
    
    findError(rowCount,columnCount-1,results, beginKick); // encontra o erro, e ja coloca na matriz!
    
    // atualizando os valores que serao utilizados na proxima iteracao!
    // so pra garantir!
    for(i=0; i < linearSystem->m_systemMatrix->m_rows;i++)
	beginKick[i] = kickBackup[i];
    
    // verificar a convergencia! e caso seja < que a tolerancia, sair do loop
     // colocando os novos valores para serem utilizados na proxima iteracao
    
     data = results->m_data[rowCount][(results->m_columns-1)];
     
//      printf("%.10f < %.10f \n", data, errMin);
     
     if ( data < errMin ){
       printf("Solucao encontrada!\n");
       break;
     }
      rowCount+=1;
   }//final do loop de iteracoes  
   
    CLEAR_SCREEN
    printf("Resultados das iteracoes:\n");
    
//     showMatrix(results);
    rowCount == iterator ? rowCount-=1, k-=1: rowCount;
    printf("Solucao: ");
    
    for(i=0; i < results->m_columns-1; i++){
      printf("X[%d] = %.10f ", i, results->m_data[rowCount][i]);
    } 
    printf("Erro: %.10f\n Total de iteracoes: %d|\n", results->m_data[rowCount][i], k+1);
    (results->m_data[rowCount][i] < errMin) ? puts("Convergiu!") : puts("Nao convergiu!");
  
  MatrixDelete(results);
  free(beginKick);
  free(kickBackup);
  
  printf("Tolerancia %.10f\n", errMin);
  return;
}

// trocar de linha vai ter de verificar se foi possivel trocar, e retornar verdadeiro se sim, ou false se nao!
static void trocarLinha(const ulong a, const ulong b, Matrix* matrix)
{
  unsigned long i;
  double *backup = (double*) malloc(matrix->m_columns* sizeof(double));
  
  //printf("Trocando |L%d = L%d|\n",a,b);
  for(i=0; i < matrix->m_columns; i++){
    backup[i] = matrix->m_data[a][i];
    matrix->m_data[a][i] = matrix->m_data[b][i];
  }
  
  for(i=0; i < matrix->m_columns; i++)
    matrix->m_data[b][i] = backup[i];
  
  //showMatrix(matrix); // only Debug
  free(backup);
  backup = NULL;
  return;
}
static int transformarEmUm(const ulong line, double k, Matrix* matrix)
{
  unsigned long j;
  if (!k)
    return FALSE;
  
 // printf("L%lu = L%lu/ %lf\n",line,line, k);
  
  for(j=0; j < matrix->m_columns; j++){
    matrix->m_data[line][j]/=k;
     if (matrix->m_data[line][j] == (-0) )
      matrix->m_data[line][j]=0;
  }
 // showMatrix(matrix); // only Debug
  return TRUE;
}
// linha onde vai comecar a zerar, linha pivo, constante, e a matriz
static int transformarEmZero(const ulong line, const ulong pivotRow, /*const int column,*/ double k, Matrix* matrix)
{
  unsigned long i;
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
   
//   puts("Colocando MAtriz e saindo");
  //showMatrix(matrix); // only Debug
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
  
  for(i=0; i < matrix->m_columns-1;i++)
    printf("Os erros encontrados foram %.10f\n",erros[i]);
  
  puts("\n\n");
  
  // encontrando o maior erro! pra armazenar na matriz!
  data = erros[0];
  for(i=1; i < matrix->m_columns-1; i++){
    if (data < erros[i]){
      data = erros[i];
      printf("O maior erro encontrado foi %.10f\n", data);
    }
  }
  matrix->m_data[line][column+1] = data;
  free(erros);
  return;
}