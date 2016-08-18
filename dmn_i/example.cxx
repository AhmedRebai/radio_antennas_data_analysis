#include "dmn_i.h"

#define N_F_DIM 10

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void myFuncF( int& N, double* X, int& NF, double* F, int* UIPARM, double* URPARM, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Function to minimise is: F( x1, ..., xN ) = sum{ ( xi - i )^2, i }
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  F[ 0 ] = 0.;
  for ( int i = 0; i < N_F_DIM; i++ ) F[ 0 ]+= ( X[ i ] - (double)i )*( X[ i ] - (double)i );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void myFuncG( int& N, double* X, int& NF, double* G, int* UIPARM, double* URPARM, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Analytical Gradient of the function to minimise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{   
  for ( int i = 0; i < N_F_DIM; i++ ) G[ i ]= 2.*( X[ i ] - (double)i );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void myFuncGH( int& N, double* X, int& NF, double* G, double* H, int* UIPARM, double* URPARM, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Analytical Gradient and Hessian of the function to minimise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  for ( int i = 0; i < N_F_DIM; i++ ) G[ i ]= 2.*( X[ i ] - (double)i );
  
  int k = 0;
  for ( int i = 0; i < N_F_DIM; i++ ) for ( int j = 0; j <= i; j++ ) {
    if ( i == j ) H[ k ] = 2.;
    else H[ k ] = 0.;
    k++; 
  }
}

void printX( std::string str, DMN_i* d )
{
  std::cout << str << ":"<< std::endl;
  for ( int i = 1; i <= N_F_DIM; i++ ) printf( "X(%2d) = %10.3lf\n", i, d->X( i ) );
}

void printX( std::string str, const std::vector<double>& v ) 
{
  std::cout << str << ":"<< std::endl;
  for ( int i = 0; i < N_F_DIM; i++ ) printf( "X(%2d) = %10.3lf\n", i+1, v[ i ] );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //                      Part I: unbounded fit
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Get a DMN toolbox
  //===
  DMN_i* xDMN = new DMN_i( N_F_DIM );
  
  // Steer DMN algorithms 'a la FORTRAN'
  //===
  xDMN->IV( 17 ) = 1000;  // Max func. evaluations
  xDMN->IV( 18 ) = 1000;  // Max iterations
  xDMN->IV( 21 ) =    0;  // Mute stdio printout
  
  // Function to minimize and its Gradient, Hessian
  //===
  xDMN->CALCF   = &myFuncF;
  xDMN->CALCG   = &myFuncG;   // Gradient
  xDMN->CALCGH  = &myFuncGH;  // Gradient and Hessian
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNF
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  // Initial Guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) xDMN->X( i ) = 10.;
  
  // Minimisation
  //===
  if ( xDMN->DMNF() ) printX( "DMNF SOLUTION ", xDMN );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNG: using analytical Gradient
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  // Initial guess as std::vector<double>
  //===
  std::vector<double> X( N_F_DIM, 10. );
  xDMN->setX( X );
  
  // Minimisation
  //===
  if ( xDMN->DMNG() ) printX( "DMNG SOLUTION ", xDMN );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNH: using analytical Hessian
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Reset initial guess
  //===
  X.clear();
  X.resize( N_F_DIM, 10. );
  
  // Minimisation using std::vector I/O
  //===
  if ( xDMN->DMNH( X ) ) printX( "DMNH SOLUTION ", X );
   
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //                      Part II: bounded fit
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  // Steer using higher level toolbox options 
  //===
  xDMN->maxFunctionEval( 100 ); // Max func. evaluations
  xDMN->maxIteration( 100 );    // Max iterations
  xDMN->printingOff();          // Mute stdio printout
  
  // Set some boundings
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) xDMN->B( 1, i ) = 2.5; // xi min value
  for( int i = 1; i <= N_F_DIM; i++ ) xDMN->B( 2, i ) = 5.5; // xi max value
   
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  // Solve with DMNFB: bounded minimisation
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  // Initial guess
  //===
  for( int i = 1; i <= N_F_DIM; i++ ) xDMN->X( i ) = 10.;
  
  // Minimisation
  //===
  if ( xDMN->DMNFB() ) printX( "DMNFB SOLUTION", xDMN );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNGB: bounded and using analytical Gradient
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Set boundings from std::vector<>
  //===
  std::vector< std::vector<double> > B( 2, std::vector<double>( N_F_DIM, 0 ) );
   
  for( int i = 0; i < N_F_DIM; i++ ) B[ 0 ][ i ] = 3.5; // xi min value
  for( int i = 0; i < N_F_DIM; i++ ) B[ 1 ][ i ] = 7.5; // xi max value
  xDMN->setB( B );
  
  // Reset initial guess
  //===
  X.clear();
  X.resize( N_F_DIM, 10. );
  
  // Minimise
  //===
  if ( xDMN->DMNGB() ) printX( "DMNGB SOLUTION", xDMN );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Solve with DMNHB: bounded and using analytical Hessian
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  // Reset initial guess
  //===
  X.clear();
  X.resize( N_F_DIM, 10. );
  
  // Minimisation using std::vector I/O
  //===
  if ( xDMN->DMNHB( X, B ) ) printX( "DMNHB SOLUTION", X );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Release the toolbox
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  delete xDMN;
  
  return 0;
}
