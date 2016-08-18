#ifndef TDOA_H
#define TDOA_H 1

#include "dmn_i/dmn_i.h"
#include <math.h>

// Maximum number of antennas
//===
#define MAX_ANTENNA 100  // Note that the memory space for input antennas data is 'static'. We allocate more space 
                         // than we require. This is convenient since the minimiser is written in FORTRAN.
			 
#define N_ANTENNA_DATA 5 // Number of antenna related data: position x,y,z, signal arrival time t and signal amplitude S.

#define ADD_USER_SPACE 3 // Additional user space for fit data. e.i. to store the cascade direction in the amplitude fit. 
			 
// Enumerations
//===
enum fitType     { POINT_SOURCE=1, PLANE_WAVE=2, EXPONENTIAL_AMPLITUDE=4 };

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			 
// Toolbox class for Point Source Fitting (TDoA)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class TDoA
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
public:
  			TDoA();
			~TDoA();

  // References to Point Source parameters
  //===
  double&		xs() { return _xs; }  // source x position
  double&		ys() { return _ys; }  // source y position
  double&		zs() { return _zs; }  // source z position
  double&		ts() { return _ts; }  // source emitting time
  
  // References to Plane Wave parameters
  //===
  double&		theta() { return _theta; }  // wave vector theta angle 
  double&		phi()   { return _phi;   }  // wave vector phi angle
  
  // Reference to wave speed
  //===
  double&		cr()    { return _cr;    }  // wave normalised speed
  
  // References to Exponential Amplitude Model
  //===
  double&		x0() { return _x0; }  // Impact x coordinate 
  double&		y0() { return _y0; }  // Impact y coordinate
  double&		z0() { return _z0; }  // Impact z coordinate: not in the fit; defined as barycenter of za values.
  double&		s0() { return _s0; }  // Signal source amplitude (in dB)
  double&		a0() { return _a0; }  // Signal characteristic loss (in dB/m: multiply by ln(10)/20 to get exponential loss factor in 1/m)
  
  // References to antennas
  //===
  int&			Na() { return _Na; }  // Number of antennas for this minimisation
  double&		xa( const int& n );   // nth antenna x position 
  double&		ya( const int& n );   // nth antenna y position
  double&		za( const int& n );   // nth antenna z position
  double&		ta( const int& n );   // nth antenna signal arrival time
  double&		sa( const int& n );   // nth antenna signal amplitude ( in dB ) 
  
  // I/O to fit
  //===
  fitType&		fitModel()   { return _fitModel;   } // Model to fit
  bool&			fixedSpeed() { return _fixedSpeed; } // Is speed fixed or is it a free parameter in the fit?
  
  // Fit methods
  //===
  double		fit( int na=-1 );    // Single fit
  double		scan( int na=-1 );   // Multiple fit with a scan of parameter space for initial conditions

private:
  DMN_i			_xDMN; // C++ interface to the FORTRAN minimiser 
  
  fitType		_fitModel;
  bool			_fixedSpeed;

  double		_Xa[ N_ANTENNA_DATA*MAX_ANTENNA + ADD_USER_SPACE ];
  int			_Na;
  
  double		_dummyReal;
  
  void			_initFit();
  void			_copyFitParameters();
  int			_nParameters;
  
  double		_scan_TDoA();
  double		_scan_EAM();
  
  double		_xs;
  double		_ys;
  double		_zs;
  double		_ts;
  double		_cs;
  
  double		_theta;
  double		_phi;
  double		_cr;
  
  double		_x0;
  double		_y0;
  double		_z0;
  double		_s0;
  double		_a0;
};

// Objective functions to minimise and their gradients
//===
void PSF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Point source
void PSF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void PWF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Plane wave
void PWF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void EAM_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Exponential amplitude model
void EAM_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );

#endif
