#include "dmn_i.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DMN_i::DMN_i( int n ):
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_dummyInt( 0 ),
_dummyReal( 0. ),
_verbosity( 1 ),
_maxFunctionEval( -1 ),
_maxIteration( -1 ),
_printingUnit( -1 ),
_LIV( 59+3*DMN_N_MAX_DIM ),
_LV( 78+DMN_N_MAX_DIM*(DMN_N_MAX_DIM+32) )
{
  reset( n );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DMN_i::~DMN_i()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DMN_i::DMN_i( const DMN_i& d )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _N = d._N;
  CALCF  = d.CALCF;
  CALCG  = d.CALCG;
  CALCGH = d.CALCGH;
  for ( int i = 0; i < _N; i++ ) {
    _X[ i ] = d._X[ i ];
    _D[ i ] = d._D[ i ];
  }
  for ( int i = 0; i < 2*_N; i++ ) _B[ i ] = d._B[ i ];
  for ( int i = 0; i < _LIV; i++ ) _IV[ i ] = d._IV[ i ];
  for ( int i = 0; i < _LV; i++ ) _V[ i ] = d._V[ i ];
  _verbosity       = d._verbosity;
  _maxFunctionEval = d._maxFunctionEval;
  _maxIteration    = d._maxIteration;
  _printingUnit    = d._printingUnit;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::reset( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _N = n;
 
  init();
  CALCF  = NULL;
  CALCG  = NULL;
  CALCGH = NULL;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::init( int alg )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define MAX_R 1e10 
{
  if ( _N == 0 ) return;  
  
  for ( unsigned  i = 0; i < _N; i++ ) { 
    _D[ i ] = 1.; 
    _X[ i ] = 0.;
    B( 1, i+1 ) = -MAX_R;
    B( 2, i+1 ) =  MAX_R; 
  }
  
  _IV[ 0 ] = 0;
  divset_( alg, _IV, _LIV, _LV, _V );
  
  // Override with user settings if required
  //===
  if ( _maxFunctionEval >= 0 ) _IV[ 16 ] = _maxFunctionEval;
  if ( _maxIteration >= 0    ) _IV[ 17 ] = _maxIteration;
  if ( _printingUnit >= 0    ) _IV[ 20 ] = _printingUnit;
}
#undef MAX_R

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& DMN_i::X( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _N ) {
    warning( "access to X is out of range" );
    _dummyReal = 0.;
    return _dummyReal;
  }
  return _X[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& DMN_i::D( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _N ) {
    warning( "access to D is out of range" );
    _dummyReal = 0.;
    return _dummyReal;
  }
  return _D[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& DMN_i::B( int i, int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 || n > _N ) || ( i <= 0 || i > 2 ) ) {
    warning( "access to B is out of range" );
    _dummyReal = 0.;
    return _dummyReal;
  }
  return _B[ 2*( n - 1 ) + i - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& DMN_i::V( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _LV ) {
    warning( "access to V is out of range" );
    _dummyReal = 0.;
    return _dummyReal;
  }
  return _V[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int& DMN_i::IV( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _LIV ) {
    warning( "access to IV is out of range" );
    _dummyInt = 0;
    return _dummyInt;
  }
  return _IV[ n - 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::maxFunctionEval( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n < 0 ) {
    warning( "Negative number of function evaluations specified." );
    return;
  }

  _maxFunctionEval = n;
  _IV[ 16 ] = n;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::maxIteration( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n < 0 ) {
    warning( "Negative number of iterations specified." );
    return;
  }

  _maxIteration = n;
  _IV[ 17 ] = n;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::printingOn()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _printingUnit = 6;
  _IV[ 20 ] = 6;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::printingOff()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _printingUnit = 0;
  _IV[ 20 ] = 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::printingUnit( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n < 0 ) {
    warning( "Negative printing unit specified." );
    return;
  }

  _printingUnit = n;
  _IV[ 20 ] = n;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DMN_i::DMNF( int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( CALCF == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( Neff == -1 ) Neff = _N;
  if ( Neff > _N ) {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }
  
  _IV[ 0 ] = 12;
  dmnf_( Neff, _D, _X, CALCF, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM );
  
  return true; 
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DMN_i::DMNG( int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( CALCF == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( CALCG == NULL ) {
    warning( "no analytical Gradient G provided. Aborting" );
    return false;
  }
  if ( Neff == -1 ) Neff = _N;
  if ( Neff > _N ) {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  _IV[ 0 ] = 12;
  dmng_( Neff, _D, _X, CALCF, CALCG, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DMN_i::DMNH( int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( CALCF == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( CALCGH == NULL ) {
    warning( "no analytical Gradient and Hessian GH provided. Aborting" );
    return false;
  }
  if ( Neff == -1 ) Neff = _N;
  if ( Neff > _N ) {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  _IV[ 0 ] = 12;
  dmnh_( Neff, _D, _X, CALCF, CALCGH, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DMN_i::DMNFB( int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( CALCF == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( Neff == -1 ) Neff = _N;
  if ( Neff > _N ) {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }
  
  _IV[ 0 ] = 12;
  dmnfb_( Neff, _D, _X, (double*)_B, CALCF, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM );
  
  return true;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DMN_i::DMNGB( int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{ 
  if ( CALCF == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( CALCG == NULL ) {
    warning( "no analytical Gradient G provided. Aborting" );
    return false;
  }
  if ( Neff == -1 ) Neff = _N;
  if ( Neff > _N ) {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  _IV[ 0 ] = 12;
  dmngb_( Neff, _D, _X, (double*)_B, CALCF, CALCG, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DMN_i::DMNHB( int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( CALCF == NULL ) {
    warning( "no F function to minimise. Aborting" );
    return false;
  }
  if ( CALCGH == NULL ) {
    warning( "no analytical Gradient and Hessian GH provided. Aborting" );
    return false;
  }
  if ( Neff == -1 ) Neff = _N;
  if ( Neff > _N ) {
    warning( "Space dimension to minimise is greater than the parameter space. Aborting" );
    return false;
  }

  _IV[ 0 ] = 12;
  dmnhb_( Neff, _D, _X, (double*)_B, CALCF, CALCGH, _IV, _LIV, _LV, _V, UIPARM, URPARM, UFPARM );
  
  return true;  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double DMN_i::valueG( int n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n <= 0 || n > _N ) {
    warning( "wrong index value while acessing to gradient G." );
    return 0.;
  }

  return( _V[ _IV[ 27 ]-1 + ( n - 1 ) ] );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::warning( const char* strwar )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( _verbosity < 1 ) return;
  
  std::cout << "DMN_i WARNING: " << strwar << std::endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::setX( const std::vector<double>& x )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( x.size() != _N ) {
    warning( "vector size does not match N while setting X. X left unchanged." );
    return;
  }
  for ( int i = 0; i < _N; i++ ) _X[ i ] = x[ i ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::setD( const std::vector<double>& d )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( d.size() != _N ) {
    warning( "vector size does not match N while setting D. D left unchanged." );
    return;
  }
  for ( int i = 0; i < _N; i++ ) _D[ i ] = d[ i ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::setB( const std::vector< std::vector< double > >& b )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( b.size() == 0 ) {
    warning( "null vector size while setting B. B left unchanged." );
    return;
  }
  
  if ( b.size() != 2 || b[ 0 ].size() != _N ) {
    warning( "vector size does not match (2,N) while setting B. B left unchanged." );
    return;
  }
  
  for ( int i = 0; i < _N; i++ ) {
    this->B( 1, i+1 ) = b[ 0 ][ i ];
    this->B( 2, i+1 ) = b[ 1 ][ i ];
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void DMN_i::getX( std::vector<double>& x )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  x.resize( _N );
  for ( int i = 0; i < _N; i++ ) x[ i ] = _X[ i ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double DMN_i::F( int* UIPARM, double* URPARM, void (*UFPARM)() )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double f0 = 0.;
  int   nf  = 0;
    
  CALCF( _N, _X, nf, &f0, UIPARM, URPARM, &UFPARM );
  
  return f0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bool DMN_i::DMNF( std::vector<double>& x, const std::vector<double>& d, int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  this->setX( x );
  if( d.size() == 0 ) {
    std::vector<double> d1( _N, 1 );
    this->setD( d1 );
  }
  else this->setD( d );
  
  if ( DMNF( UIPARM, URPARM, UFPARM, Neff ) ) {
    this->getX( x );
    return true;
  }
  return false;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bool DMN_i::DMNG( std::vector<double>& x, const std::vector<double>& d, int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  this->setX( x );
  if( d.size() == 0 ) {
    std::vector<double> d1( _N, 1 );
    this->setD( d1 );
  }
  else this->setD( d );
  
  if ( DMNG( UIPARM, URPARM, UFPARM, Neff ) ) {
    this->getX( x );
    return true;
  } 
  return false;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bool DMN_i::DMNH( std::vector<double>& x, const std::vector<double>& d, int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  this->setX( x );
  if( d.size() == 0 ) {
    std::vector<double> d1( _N, 1 );
    this->setD( d1 );
  }
  else this->setD( d );
  
  if ( DMNH( UIPARM, URPARM, UFPARM, Neff ) ) {
    this->getX( x );
    return true;
  }
  return false;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bool DMN_i::DMNFB( std::vector<double>& x, const std::vector< std::vector<double> >& b, const std::vector<double>& d, int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  this->setX( x );
  this->setB( b );
  if( d.size() == 0 ) {
    std::vector<double> d1( _N, 1 );
    this->setD( d1 );
  }
  else this->setD( d );
  
  if ( DMNFB( UIPARM, URPARM, UFPARM, Neff ) ) {
    this->getX( x );
    return true;
  }  
  return false;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bool DMN_i::DMNGB( std::vector<double>& x, const std::vector< std::vector<double> >& b, const std::vector<double>& d, int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  this->setX( x );
  this->setB( b );
  if( d.size() == 0 ) {
    std::vector<double> d1( _N, 1 );
    this->setD( d1 );
  }
  else this->setD( d );
  
  if ( DMNGB( UIPARM, URPARM, UFPARM, Neff ) ) {
    this->getX( x );
    return true;
  } 
  return false;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
bool DMN_i::DMNHB( std::vector<double>& x, const std::vector< std::vector<double> >& b, const std::vector<double>& d, int* UIPARM, double* URPARM, void (*UFPARM)(), int Neff )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  this->setX( x );
  this->setB( b );
  if( d.size() == 0 ) {
    std::vector<double> d1( _N, 1 );
    this->setD( d1 );
  }
  else this->setD( d );
  
  if ( DMNHB( UIPARM, URPARM, UFPARM, Neff ) ) {
    this->getX( x );
    return true;
  } 
  return false;
}
