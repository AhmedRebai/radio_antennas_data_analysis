#include "TDoA.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TDoA::TDoA():
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Class constructor. Main task is tuning the DMN_i interfaces handling 
//  the minimisation.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_xDMN( 5 ) // There are up to 5 parameters to minimise: the source position (3) + the emission time (1) + the wave speed.
{
  // Initialisation of antennas times and positions
  //===
  for ( int i = 0; i < N_ANTENNA_DATA*MAX_ANTENNA + ADD_USER_SPACE; i++ ) _Xa[ i ] = 0.0;
  _Na = 0;
  
  // Minimiser tuning
  //===
  _xDMN.maxFunctionEval( 1000 );
  _xDMN.maxIteration( 1000 );
  _xDMN.printingOff(); // Mute the minimiser, comment this for debug
  
  // Default fit settings
  //===
  _fitModel    = POINT_SOURCE;
  _fixedSpeed  = true;
  _cr          = 1.0;
  
  // Default parameter initialisation
  //===
  _xs = _ys = _zs = _ts = 0.0;
  _theta = _phi = 0.0;
  _x0 = _y0 = _z0 = _s0 = _a0 = 0.0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TDoA::~TDoA()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TDoA::_initFit()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define MAX_R 1e10 
{
  // Reset constrains
  //===
  for ( int i = 0; i < _xDMN.N(); i++ ) { 
    _xDMN.B( 1, i+1 ) = -MAX_R;
    _xDMN.B( 2, i+1 ) =  MAX_R; 
  }

  if ( _fitModel == POINT_SOURCE ) {
  
    _nParameters = 4;

    // Connect the objective function to minimise and its gradient 
    //===
    _xDMN.CALCF = &PSF_function_F;
    _xDMN.CALCG = &PSF_function_G;
  
    // Set bounding: Ts <= 0
    //===
    _xDMN.B( 2, 4 ) = 0.0; // Ts parameter (4th) max value (2) is 0.0
    
    // Initialise fit parameters
    //===
    _xDMN.X( 1 ) = _xs;
    _xDMN.X( 2 ) = _ys;
    _xDMN.X( 3 ) = _zs;
    _xDMN.X( 4 ) = _ts;
    _xDMN.X( 5 ) = _cr; 
  }
  else if ( _fitModel == PLANE_WAVE ) {
  
    _nParameters = 2;

    // Connect the objective function to minimise and its gradient 
    //===
    _xDMN.CALCF = &PWF_function_F;
    _xDMN.CALCG = &PWF_function_G;
    
    // Set bounding: 0 <= theta <= pi
    //===
    _xDMN.B( 1, 1 ) = 0.0;    // theta parameter (1st) min value (1) is 0.0
    _xDMN.B( 2, 1 ) = M_PI;   // theta parameter (1st) max value (2) is pi
    _xDMN.B( 1, 2 ) = 0.0;    // phi parameter (2nd) min value (1) is 0.0
    _xDMN.B( 2, 2 ) = 2*M_PI; // phi parameter (2nd) max value (2) is 2*pi
    
    // Initialise fit parameters
    //===
    _xDMN.X( 1 ) = _theta;
    _xDMN.X( 2 ) = _phi;
    _xDMN.X( 3 ) = _cr;
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE ) {
  
    _nParameters = 4;

    // Connect the objective function to minimise and its gradient 
    //===
    _xDMN.CALCF = &EAM_function_F;
    _xDMN.CALCG = &EAM_function_G;
    
    // Set bounding: a0 >= 0
    //===
    _xDMN.B( 1, 4 ) = 0.0;    // a0 parameter (4th) min value (1) is 0.0
    
    // Initialise fit parameters
    //===
    _xDMN.X( 1 ) = _x0;
    _xDMN.X( 2 ) = _y0;
    _xDMN.X( 3 ) = _s0;
    _xDMN.X( 4 ) = _a0;
    
    // Initialise cascade direction
    //===
    _Xa[ N_ANTENNA_DATA*MAX_ANTENNA     ] = cos( _phi )*sin( _theta );
    _Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 1 ] = sin( _phi )*sin( _theta );
    _Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 2 ] = cos( _theta );
  }
  
  if ( !_fixedSpeed && ( _fitModel == POINT_SOURCE || _fitModel == PLANE_WAVE ) ) {
    
    //  Add wave speed as last parameter
    //===
    _nParameters++;
    
    // Set bounding: cr >= 0
    //===
    _xDMN.B( 1, _nParameters ) = 0.0; // Speed parameter min value is 0.0
  }
}
#undef MAX_R

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TDoA::_copyFitParameters()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
{
  if ( _fitModel == POINT_SOURCE ) {
    // Get fit parameters
    //===
    _xs = _xDMN.X( 1 );
    _ys = _xDMN.X( 2 );
    _zs = _xDMN.X( 3 );
    _ts = _xDMN.X( 4 );
    if ( _nParameters == 5 ) _cr = _xDMN.X( 5 );
  }
  else if ( _fitModel == PLANE_WAVE ) {
    // Get fit parameters
    //===
    _theta = _xDMN.X( 1 );
    _phi   = _xDMN.X( 2 );
    _cr    = _xDMN.X( 3 );
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE ) {
    // Get fit parameters
    //===
    _x0 = _xDMN.X( 1 );
    _y0 = _xDMN.X( 2 );
    _s0 = _xDMN.X( 3 );
    _a0 = _xDMN.X( 4 );
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& TDoA::xa( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna x position.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n > MAX_ANTENNA ) {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else return _Xa[ N_ANTENNA_DATA*( n - 1 )    ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& TDoA::ya( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna y position.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n > MAX_ANTENNA ) {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 1 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& TDoA::za( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna z position.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n > MAX_ANTENNA ) {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 2 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& TDoA::ta( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna signal arrival time.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n > MAX_ANTENNA ) {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 3 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& TDoA::sa( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna signal amplitude (in dB).
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( n > MAX_ANTENNA ) {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 4 ];
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double TDoA::fit( int na )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Perform a single fit with user suplied initial conditions. Returns the 
//  fit chi2. In case of faillure -1 is returned. The fit best guess is 
//  stored in fit parameters vector.
//
//  The number na of antenna can be provided as input argument. The current
//  value of _Na is assumed otherwise. 
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  // Update number of antennas if required
  //===
  if ( na > 0 ) _Na = na;
  
  double t0, z0;
  if ( _fitModel == POINT_SOURCE || _fitModel == PLANE_WAVE ) {
    // Regularise the time origin
    //===
    t0 = this->ta( 1 );
    for ( int i = 2; i <= _Na; i++ ) if ( this->ta( i ) < t0 ) t0 = this->ta( i );
    for ( int i = 1; i <= _Na; i++ ) this->ta( i )-= t0;
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE ) {
    // Regularise z coordinate
    //===
    z0 = 0.0;
    for ( int i = 1; i <= _Na; i++ ) z0+= this->za( i )/_Na;
    for ( int i = 1; i <= _Na; i++ ) this->za( i )-= z0;
    this->z0() = z0;
  }
  
  // Initialise according to fit model
  //===
  _initFit();
  
  // Do the fit
  //===
  double chi2 = -1;
  if ( _xDMN.DMNGB( &_Na, _Xa, NULL, _nParameters ) ) {
    chi2 = _xDMN.valueF(); 
    this->_copyFitParameters();
  }
  
  if ( _fitModel == POINT_SOURCE || _fitModel == PLANE_WAVE ) {
    // Restore the time origin
    //===
    for ( int i = 1; i <= _Na; i++ ) this->ta( i )+= t0;
    if ( _fitModel == POINT_SOURCE ) {
      this->ts()+= t0;
      _xDMN.X( 4 )+= t0;
    }
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE ) {
    // Restore the z origin
    //===
    for ( int i = 1; i <= _Na; i++ ) this->za( i )+= z0;
  }
  
  return chi2;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double TDoA::scan( int na )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Generic routing for scan like algorithms. Scans try several fit with various
//  initial conditions. The minimum chi2 over all fits is returned. In case of 
//  faillure -1 is returned. The fit best guess is stored as parameter vector.
//
//  The number na of antenna can be provided as input argument. The current
//  value of _Na is assumed otherwise.  
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  // Update number of antennas if required
  //===
  if ( na > 0 ) _Na = na;
  
  // Initialise according to fit model
  //===
  _initFit();

  // Route on specific scan algorithm
  //===
  if ( _fitModel == POINT_SOURCE || _fitModel == PLANE_WAVE ) {
    return this->_scan_TDoA();
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE ) {
    return this->_scan_EAM(); 
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double TDoA::_scan_TDoA()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Scan Algorithm for TDoA algorithms. Tries several fit with various initial
//  conditions.  For example for point source localisation initial conditions 
//  are taken as a scan of the space in spherical coordinates: r, theta, phi. 
//  Steps by 15 degrees in theta, half-quadrant (45 degrees) in phi and 
//  logarithmically in r, from 1 m to 10 km.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define N_SPEED  5
#define N_RANGE  5
#define N_THETA  8
#define N_PHI    9
{ 
  // Scan values
  //===
  double speed_v[] = {  0.20, 0.40, 0.60, 0.80, 1.0  };
  double range_v[] = {  1., 10., 100., 1000., 10000. };
  double theta_v[] = {  -15., 0., 15., 30., 45., 60., 75., 90. };
  double phi_v[]   = {  0.,  45.,  90., 135., 180., 225., 270., 315., 360. };
  
  double deg2rad = M_PI/180.;
  
  // Tune scan parameters according to wave type fit
  //===
  double x0, y0, z0;
  int ir_max, ic_max;
  if ( _fitModel == POINT_SOURCE ) {
    x0 = y0 = z0 = 0.0;
    for ( int i = 1; i <= _Na; i++ ) {
      x0+= this->xa( i )/_Na;
      y0+= this->ya( i )/_Na;
      z0+= this->za( i )/_Na;
    }
    ir_max = N_RANGE;
  }
  else if ( _fitModel == PLANE_WAVE ) {
    ir_max = 1;
  }
  
  if ( _fixedSpeed ) ic_max = 1;
  else ic_max = N_SPEED;
  
  // Do the scan
  //==
  double chi2 = -1;
  std::vector<double> S( _nParameters, 0 );
  for ( int ic = 0; ic < ic_max; ic++ ) for ( int ir = 0; ir < ir_max; ir++ ) for ( int it = 0; it < N_THETA; it++ ) for ( int ip = 0; ip < N_PHI; ip++ ) {
  
    if ( !_fixedSpeed ) this->cr() = speed_v[ ic ];
  
    if ( _fitModel == POINT_SOURCE ) {
      this->xs() = range_v[ ir ]*cos( phi_v[ ip ]*deg2rad )*sin( theta_v[ it ]*deg2rad ) + x0;
      this->ys() = range_v[ ir ]*sin( phi_v[ ip ]*deg2rad )*sin( theta_v[ it ]*deg2rad ) + y0;
      this->zs() = range_v[ ir ]*cos( theta_v[ it ]*deg2rad ) + z0;
      this->ts() = -range_v[ ir ]/this->cr();
    }
    else if ( _fitModel == PLANE_WAVE ) {
      this->theta() = theta_v[ it ]*deg2rad;
      this->phi()   = phi_v[ ip ]*deg2rad;
    }
    
    double d = this->fit();
    if ( ( d >= 0 && d < chi2 ) || ( chi2 < 0 ) ) {
      chi2   = d;
      for ( int ipar = 0; ipar < _nParameters; ipar++ ) {
        S[ ipar ] = _xDMN.X( ipar+1 );  
      }
    }  
  }
  
  // Restore best parameter guess
  //===
  for ( int ipar = 0; ipar < _nParameters; ipar++ ) {
    _xDMN.X( ipar+1 ) = S[ ipar ];
  }
  this->_copyFitParameters(); 
  
  return chi2;
}
#undef N_PHI
#undef N_THETA
#undef N_RANGE
#undef N_SPEED

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double TDoA::_scan_EAM()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Scan algorithm for amplitude reconstruction model according to exponential
//  loss with the lateral distance to the shower axis. Initial conditions are
//  taken as a scan of the impact point over the field of antennas. Loss is
//  scanned for a characteristic length in the range 50-500m.    
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define N_RANGE  11
#define N_LAMBDA 7
{
  // Get array half widths along x and y directions and get the coordinates of 
  // the closest antenna to shower impact point 
  //===
  double dmin, dmax, Dx, Dy, x0, y0, s0;
  int i0;
  
  // Half width along x
  //===
  dmin = dmax =  this->xa( 1 ); 
  for ( int i = 2; i <= _Na; i++ ) {
    if ( this->xa( i ) < dmin ) dmin = this->xa( i );
    else if ( this->xa( i ) > dmax ) dmax = this->xa( i );
  }
  Dx = 0.5*( dmax - dmin );
  
  // Half width along y 
  //===
  dmin = dmax =  this->ya( 1 );
  for ( int i = 2; i <= _Na; i++ ) {
    if ( this->ya( i ) < dmin ) dmin = this->ya( i );
    else if ( this->ya( i ) > dmax ) dmax = this->ya( i );
  }
  Dy = 0.5*( dmax - dmin ); 
  
  // Coordinates of closest antenna to shower impact point
  //===
  i0 = 1;
  dmax =  this->sa( 1 );
  for ( int i = 2; i <= _Na; i++ ) {
    if ( this->sa( i ) > dmax ) {
      dmax = this->sa( i );
      i0 = i;
    }
  }
  x0 = this->xa( i0 );
  y0 = this->ya( i0 );
  s0 = this->sa( i0 );
    
  // Scan values
  //===
  double lambda_v[] = {  50., 100., 150., 200., 300., 400., 500. };
  double range_v[]  = {  -4.0, -2.0, -1.0, -0.50, -0.25, 0.0, 0.25, 0.50, 1.0, 2.0, 4.0 };
  
  // Do the scan
  //===
  double chi2 = -1;
  std::vector<double> S( _nParameters, 0 );
  for ( int il = 0; il < N_LAMBDA; il++ ) for ( int ix = 0; ix < N_RANGE; ix++ ) for ( int iy = 0; iy < N_RANGE; iy++ ) {
       
    double dx = range_v[ ix ]*Dx;
    double dy = range_v[ iy ]*Dy;
    this->x0() = x0 + dx;
    this->y0() = y0 + dy;
    this->s0() = s0*exp( sqrt( dx*dx + dy*dy )/lambda_v[ il ] );
    this->a0() = 8.6859/lambda_v[ il ];
   
    double d = this->fit();
    if ( ( d >= 0 && d < chi2 ) || ( chi2 < 0 ) ) {
      chi2   = d;
      for ( int ipar = 0; ipar < _nParameters; ipar++ ) {
        S[ ipar ] = _xDMN.X( ipar+1 );
      }
    }  
  }
  
  // Restore best parameter guess
  //===
  for ( int ipar = 0; ipar < _nParameters; ipar++ ) {
    _xDMN.X( ipar+1 ) = S[ ipar ];  
  }
  this->_copyFitParameters();
  
  return chi2;
}
#undef N_RANGE
#undef N_LAMBDA

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for point source fit. It is built 
//  as a chi2 over all antennas, requiring:
//
//  | Xa - Xs |^2 - cs^2*( Ta - Ts )^2 = 0
//
//  Time is assumed to be expressed in m and speed cs is normalised to C0  
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double d, fi;
  
  F[ 0 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ ) {
    d  = ( X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ] );
    fi = d*d;
    
    d  = ( X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] );
    fi+= d*d;
    
    d  = ( X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] );
    fi+= d*d;
    
    d  = X[ 4 ]*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] );
    fi-= d*d;
    
    F[ 0 ]+= fi*fi;
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double d, fi;
  
  G[ 0 ] = G[ 1 ] = G[ 2 ] = G[ 3 ] = G[ 4 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ ) {
    d  = ( X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ] );
    fi = d*d;
    
    d  = ( X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] );
    fi+= d*d;
    
    d  = ( X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] );
    fi+= d*d;
    
    d  = X[ 4 ]*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] );
    fi-= d*d;
    
    G[ 0 ]+= ( X[ 0 ] - Xa[ i*N_ANTENNA_DATA	] )*fi;
    G[ 1 ]+= ( X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] )*fi;
    G[ 2 ]+= ( X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] )*fi;
    G[ 3 ]-= X[ 4 ]*X[ 4 ]*( X[ 3 ] - Xa[ N_ANTENNA_DATA + 3 ] )*fi;
    if ( N == 5 ) G[ 4 ]-= X[ 4 ]*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] )*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] )*fi; 
  }
  G[ 0 ]*= 4;
  G[ 1 ]*= 4;
  G[ 2 ]*= 4;
  G[ 3 ]*= 4;
  if ( N == 5 ) G[ 4 ]*= 4;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PWF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for plane wave fit. It is built 
//  as a chi2 over all antennas, requiring:
//
//  ( Xa_j - Xa_i ).k - cr*( Ta_j - Ta_i ) = 0
//
//  over all pairs.
//
//  Time is assumed to be expressed in m and speed cr is normalised to C0  
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double st, ct, sp, cp, fi;
  
  ct = cos( X[ 0 ] );
  st = sin( X[ 0 ] );
  cp = cos( X[ 1 ] );
  sp = sin( X[ 1 ] );
  
  F[ 0 ] = 0.0;
  for ( int j = 0; j < Na[ 0 ]; j++ ) for ( int i = j+1; i < Na[ 0 ]; i++ ) {
    fi  = ( Xa[ j*N_ANTENNA_DATA     ] - Xa[ i*N_ANTENNA_DATA     ] )*st*cp;
    fi += ( Xa[ j*N_ANTENNA_DATA + 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] )*st*sp;
    fi += ( Xa[ j*N_ANTENNA_DATA + 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] )*ct;
    fi -= X[ 2 ]*( Xa[ j*N_ANTENNA_DATA + 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] );
    F[ 0 ]+= fi*fi;
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PWF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double st, ct, sp, cp, xij, yij, zij, tij, fi;
  
  ct = cos( X[ 0 ] );
  st = sin( X[ 0 ] );
  cp = cos( X[ 1 ] );
  sp = sin( X[ 1 ] );
  
  G[ 0 ] = G[ 1 ] = G[ 2 ] = 0.0;
  for ( int j = 0; j < Na[ 0 ]; j++ ) for ( int i = j+1; i < Na[ 0 ]; i++ ) {
  
    xij = Xa[ j*N_ANTENNA_DATA     ] - Xa[ i*N_ANTENNA_DATA     ];
    yij = Xa[ j*N_ANTENNA_DATA + 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zij = Xa[ j*N_ANTENNA_DATA + 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ];
    tij = Xa[ j*N_ANTENNA_DATA + 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ];
  
    fi  = ( xij*cp + yij*sp )*st + zij*ct - X[2]*tij;
    
    G[ 0 ]+= ( ( xij*cp + yij*sp )*ct - zij*st )*fi;
    G[ 1 ]+= ( ( yij*cp - xij*sp )*st + zij*ct )*fi;
    G[ 2 ]-= tij*fi;
  }
  G[ 0 ]*= 2;
  G[ 1 ]*= 2;
  G[ 2 ]*= 2;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EAM_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for exponential amplitude fit. It is built 
//  as a chi2 over all antennas, requiring:
//
//  si = s0 - a0*di
//
//  where si is the signal on antenna i and di the lateral distance to the axis
//  given as:
//
//  di^2 = | Xi - X0 |^2 - ( ( Xi - X0 ).u )^2
//
//  with X0 the point on the cascade axis that intercepts the plane z0 = 0, and
//  u the cascade direction  
//
//  Amplitude si and loss a0 are assumed to be expressed in dB and dB/m.  
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double fi, di, xi, yi, zi, ux, uy, uz;
  
  ux = Xa[ N_ANTENNA_DATA*MAX_ANTENNA     ];
  uy = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 1 ];
  uz = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 2 ];
 
  F[ 0 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ ) {
    xi = X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ];
    yi = X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zi =        - Xa[ i*N_ANTENNA_DATA + 2 ];
    di = xi*ux + yi*uy + zi*uz; 
    di = sqrt( xi*xi + yi*yi + zi*zi - di*di );
    fi = X[ 3 ]*di + Xa[ i*N_ANTENNA_DATA + 4 ] - X[ 2 ];
    F[ 0 ]+= fi*fi;
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EAM_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double fi, si, di, xi, yi, zi, ux, uy, uz;
  
  ux = Xa[ N_ANTENNA_DATA*MAX_ANTENNA     ];
  uy = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 1 ];
  uz = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 2 ];
  
  G[ 0 ] = G[ 1 ] = G[ 2 ] = G[ 3 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ ) {  
    xi = X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ];
    yi = X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zi =        - Xa[ i*N_ANTENNA_DATA + 2 ];
    si = xi*ux + yi*uy + zi*uz; 
    di = sqrt( xi*xi + yi*yi + zi*zi - si*si );
    fi = X[ 3 ]*di + Xa[ i*N_ANTENNA_DATA + 4 ] - X[ 2 ];
    
    G[ 0 ]+= fi/di*( xi - si*ux );
    G[ 1 ]+= fi/di*( yi - si*uy );
    G[ 2 ]+= fi;
    G[ 3 ]+= fi*di;
  }
  G[ 0 ]*= 2*X[ 3 ];
  G[ 1 ]*= 2*X[ 3 ];
  G[ 2 ]*= -2;
  G[ 3 ]*= 2;
}
