//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//  C++ interface to FORTRAN DMN functions
//  Required FORTRAN file: dmnfhg.f
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef DMN_I_H
#define DMN_I_H 1

#include <iostream>
#include <vector>

#define DMN_N_MAX_DIM 200

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extern "C" { 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Prototyping of FORTRAN functions
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void dmnf_(  int& N, double* D, double* X, 
               void (*CALCF)(int&, double*, int&, double*, int*, double*, void*),
	       int* IV, int& LIV, int& LV, double* V, int* UIPARM, double* URPARM, 
	       void (*UFPARM)() );
 
  void dmng_(  int& N, double* D, double* X, 
               void (*CALCF)(int&, double*, int&, double*, int*, double*, void*), 
	       void (*CALCG)(int&, double*, int&, double*, int*, double*, void*), 
	       int* IV, int& LIV, int& LV, double* V, int* UIPARM, double* URPARM, 
	       void (*UFPARM)() );
	      
  void dmnh_(  int& N, double* D, double* X, 
               void (*CALCF)(int&, double*, int&, double*, int*, double*, void*), 
	       void (*CALCGH)(int&, double*, int&, double*, double*, int*, double*, void*), 
	       int* IV, int& LIV, int& LV, double* V, int* UIPARM, double* URPARM, 
	       void (*UFPARM)() );
	      
  void dmnfb_( int& N, double* D, double* X, double* B,
               void (*CALCF)(int&, double*, int&, double*, int*, double*, void*),
	       int* IV, int& LIV, int& LV, double* V, int* UIPARM, double* URPARM, 
	       void (*UFPARM)() );
 
  void dmngb_( int& N, double* D, double* X, double* B, 
               void (*CALCF)(int&, double*, int&, double*, int*, double*, void*), 
	       void (*CALCG)(int&, double*, int&, double*, int*, double*, void*), 
	       int* IV, int& LIV, int& LV, double* V, int* UIPARM, double* URPARM, 
	       void (*UFPARM)() );
	      
  void dmnhb_( int& N, double* D, double* X, double* B, 
               void (*CALCF)(int&, double*, int&, double*, int*, double*, void*), 
	       void (*CALCGH)(int&, double*, int&, double*, double*, int*, double*, void*), 
	       int* IV, int& LIV, int& LV, double* V, int* UIPARM, double* URPARM, 
	       void (*UFPARM)() );
	      
  void divset_( int& ALG, int* IV, int& LIV, int& LV, double* V );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class DMN_i {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Interface class to DMN toolkit
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
public:

  // Constructor/Destructor
  //===
  			DMN_i( int n=0 );
  			~DMN_i();

  // Pointers to function F to minimise, Gradient G and Hessian H
  //===			
  void 			(*CALCF)( int&, double*, int&, double*, int*, double*, void*);
  void 			(*CALCG)( int&, double*, int&, double*, int*, double*, void*);
  void 			(*CALCGH)(int&, double*, int&, double*, double*, int*, double*, void*);

  // Public initialisation
  //===
  void			reset( int n=0 );
  void			init( int algorithm=2 );			

  //  Accessors to data members
  //===
  int			N()   { return _N;   } // set via init only
  int			LIV() { return _LIV; } //
  int			LV()  { return _LV;  } //

  double&		X(  int n );
  double&		D(  int n );
  double&		B(  int i, int n );
  int&			IV( int n );
  double&		V(  int n );

  // Minimisation algorithms: DMNF, DMNG and DMNH and their bounded versions
  //===
  bool			DMNF(  int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );    
  bool			DMNG(  int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  bool			DMNH(  int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  bool			DMNFB( int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );    
  bool			DMNGB( int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  bool			DMNHB( int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  
  //  Minimiser I/O
  //===
  double		F( int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL );
  int			minimisationStatus() { return _IV[ 0 ]; }
  double		valueF() { return _V[ 9 ]; }
  double		valueG( int n );
  int			nIteration()   { return _IV[ 30 ]; }
  int			nEvaluationF() { return _IV[  5 ]; }
  int			nEvaluationG() { return _IV[ 29 ]; }
  			  
  // std::vector<> interface
  //===
  void			setX(  const std::vector<double>& );
  void			setD(  const std::vector<double>& );
  void			setB(  const std::vector< std::vector< double > >& );
  void   		getX( std::vector<double>& );
  
  bool   		DMNF(  std::vector<double>& X, const std::vector<double>& D=std::vector<double>( 0 ), int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );    
  bool   		DMNG(  std::vector<double>& X, const std::vector<double>& D=std::vector<double>( 0 ), int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  bool   		DMNH(  std::vector<double>& X, const std::vector<double>& D=std::vector<double>( 0 ), int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  bool   		DMNFB( std::vector<double>& X, const std::vector< std::vector<double> >& B0, const std::vector<double>& D=std::vector<double>( 0 ), int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );    
  bool   		DMNGB( std::vector<double>& X, const std::vector< std::vector<double> >& B0, const std::vector<double>& D=std::vector<double>( 0 ), int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  bool   		DMNHB( std::vector<double>& X, const std::vector< std::vector<double> >& B0, const std::vector<double>& D=std::vector<double>( 0 ), int* UIPARM=NULL, double* URPARM=NULL, void (*UFPARM)()=NULL, int Neff=-1 );
  
  // Higher level control
  //===
  int			maxFunctionEval() { return _maxFunctionEval; }
  int			maxIteration()    { return _maxIteration;    }
  int			printingUnit()    { return _printingUnit;    }
  void			maxFunctionEval( int );
  void			maxIteration( int );
  void			printingOn();
  void			printingOff();
  void			printingUnit( int );

  
  // Warnings
  //===
  void			warningsOn()  { _verbosity = 1; }
  void			warningsOff() { _verbosity = 0; }
  void 			warning( const char* );	

private:
  int			_N;
  double		_D[ DMN_N_MAX_DIM ];
  double		_X[ DMN_N_MAX_DIM ];
  double		_B[ 2*DMN_N_MAX_DIM ];
  int			_IV[ 59+3*DMN_N_MAX_DIM ];
  int			_LIV;
  int			_LV;
  double		_V[ 78+DMN_N_MAX_DIM*(DMN_N_MAX_DIM+32) ];
  
  int			_dummyInt;
  double   		_dummyReal;
  
  char			_verbosity;
  
  int			_maxFunctionEval;
  int			_maxIteration;
  int			_printingUnit;
  
  			DMN_i( const DMN_i& ); // copy is protected
			
};

#endif
