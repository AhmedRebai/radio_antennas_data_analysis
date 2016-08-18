#include "TDoA.h"

#define FIT_TYPE 1          // 1: PLANE WAVE otherwise POINT SOURCE fit.

#define RAD2DEG 180./M_PI

int main()
{
  TDoA xTDoA; // the TDoA toolbox
  
  // Data files (I/Os)
  //===
  #if( FIT_TYPE == 1 )
    FILE* fid_in  = fopen( "planwave-data.txt",   "r" );
    FILE* fid_out = fopen( "EAS-parameters.txt", "w"  );
  #else
    FILE* fid_in  = fopen( "pointsource-data.txt", "r" );
    FILE* fid_out = fopen( "source-position.txt",  "w" );
  #endif  
  
  // Loop on antenna data
  //===
  for ( int iToy = 0; iToy < 10; iToy ++ )   // Number of events is hard-coded
  {
    // Read this event
    //===
    fscanf( fid_in, "%d", &xTDoA.Na() );     // Get/set the number of antennas in the event
    for ( int i = 1; i <= xTDoA.Na(); i++ ) 
    {
      fscanf( fid_in,                        // Get/set antennas positions, signal arrival times and level, in dB.
              "%lf%lf%lf%lf%lf",     
              &xTDoA.xa( i ), 
	      &xTDoA.ya( i ), 
	      &xTDoA.za( i ), 
	      &xTDoA.ta( i ), 
	      &xTDoA.sa( i ) );   
    } 

    // TDoA fit settings
    //===
    #if( FIT_TYPE == 1 )
      xTDoA.fitModel()  = PLANE_WAVE;   // Set either to PLANE_WAVE or POINT_SOURCE
    #else
      xTDoA.fitModel()  = POINT_SOURCE;
    #endif
    xTDoA.fixedSpeed()  = true;        // True: wave speed is a fixed parameter in the fit / false: wave speed is an additional free parameter in the fit
    xTDoA.cr()          = 0.9;         // Speed value to use if fixed. By default it is set to 1.0 at initialisation.
    
    // Do the TDoA fit
    //===
    double chi2 = xTDoA.scan();  // Perform a scan for the xTDoA.Na() first antennas
  
    // Output TDoA result
    //===
    #if( FIT_TYPE == 1 )
      fprintf( fid_out, 
               "%15.5le %15.5le %15.5le %15.5le   |", // plane wave direction 
	       xTDoA.theta()*RAD2DEG, 
	       xTDoA.phi()*RAD2DEG, 
	       xTDoA.cr(), chi2 ); 
    #else   
      fprintf( fid_out, 
               "%15.5le %15.5le %15.5le %15.5le %15.5le %15.5le\n", // point source parameters
	       xTDoA.xs(), 
	       xTDoA.ys(), 
	       xTDoA.zs(), 
	       xTDoA.ts(), 
	       xTDoA.cr(), 
	       chi2 );
    #endif
    
    #if( FIT_TYPE == 1 ) 
      // EAM fit settings
      //===
      xTDoA.fitModel()  = EXPONENTIAL_AMPLITUDE;
    
      // Do the EAM fit
      //===
      chi2 = xTDoA.scan();  // Perform a scan for the xTDoA.Na() first antennas
    
      fprintf( fid_out, 
               "%15.5le %15.5le %15.5le %15.5le %15.5le %15.5le\n", // EAS source level, lateral attenuation and impact point. 
	       xTDoA.x0(), 
	       xTDoA.y0(), 
	       xTDoA.z0(),      
               xTDoA.s0(), 
	       8.6859/xTDoA.a0(), // Convert attenuation in dB/m to attenuation length in m => factor 8.6859 = 20*log_{10}(e)
	       chi2 );  
    #endif                                                                                                                           
  }
  
  // Close I/Os
  //===
  fclose( fid_out );
  fclose( fid_in  );

  return 0;
}
