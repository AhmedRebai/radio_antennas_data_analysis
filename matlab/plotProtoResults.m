function main

  clear all;
  close all;
  
  %% Number of simulated events
  %%===
  Nevt = 1000;
  
  %% Truth
  %%===
  thetaTrue  = 20;
  phiTrue    = 60;
  X0True     = [ -4.72667e+001 1.88333e+000 2.66822e+003 ];
  S0True     = 20*log10( 1024 );
  lambdaTrue = 100;
  
  %% Read fit inputs
  %%===
  Na   = 6;
  fid  = fopen( sprintf( 'planwave-data-%03d.txt', lambdaTrue ), 'r' );
  data = fscanf( fid, '%f', [ 1+5*Na, Nevt ] );
  fclose( fid );
  Ka = [];
  for k = 0:Na-1; Ka = [ Ka; 5*k + [ 1:3 ] + 1 ]; end;
  Xa      = zeros( Na, 3 ); 
  Xa( : ) = data( Ka, 1 );
  Ta      = data( [ 5 + 5*[ 0:Na-1 ] ], : );
  Sa      = data( [ 6 + 5*[ 0:Na-1 ] ], : );

  %% Read fit results
  %%===
  fid = fopen( sprintf( 'EAS-parameters-%03d.txt', lambdaTrue ), 'r' );
  for k = 1:Nevt
    data = fscanf( fid, '%f', [ 1, 4 ] );
    theta( k )   = data( 1 );
    phi( k )     = data( 2 );
    cr( k )      = data( 3 );
    chi2PWF( k ) = data( 4 );
    data = fscanf( fid, '%s', 1 );
    data = fscanf( fid, '%f', [ 1, 6 ] );
    x0( k )      = data( 1 );
    y0( k )      = data( 2 );
    z0( k )      = data( 3 );
    s0( k )      = data( 4 );
    lambda( k )  = data( 5 );
    chi2EAM( k ) = data( 6 );
  end
  fclose( fid );
  
  %% Plot resolutions: wave direction & EAS amplitude, lateral loss
  %%===
  figure; 
  subplot( 2, 2, 1 );
  myHist( theta - thetaTrue, 20, 'k', '\Delta\theta ( deg )' );
  subplot( 2, 2, 2 );
  myHist( phi - phiTrue, 20, 'k', '\Delta\phi ( deg )' );
  subplot( 2, 2, 3 );
  X = s0 - S0True;
  K = find( abs( X ) < 20 );
  myHist( X( K ), 20, 'k', '\DeltaS0 ( dB )' );
  subplot( 2, 2, 4 );
  X = 20*log10( lambda/lambdaTrue );
  K = find( abs( X ) < 20 );
  myHist( X( K ), 20, 'k', '\Delta\lambda ( dB )' );
  
  %% Plot resolutions: EAS impact point
  %%===
  figure; 
  subplot( 1, 2, 1 );
  X = x0 - X0True( 1 );
  K = find( abs( X ) < 100 );
  myHist( X( K ), 20, 'k', '\DeltaX0 ( m )' );
  subplot( 1, 2, 2 );
  X = y0 - X0True( 2 );
  K = find( abs( X ) < 100 );
  myHist( X( K ), 20, 'k', '\DeltaY0 ( m )' );
  
  %% Plot EAS fit result ( selected event )
  %%===
  figure; 
  kevt = 1; %% selected event
  
  subplot( 1, 2, 1 ); %% Amplitude as lateral distance
  dX = [ Xa - ones( Na, 1 )*[ x0( kevt ), y0( kevt ), z0( kevt ) ] ];
  cp = cos( phi( kevt )*deg2rad(1) ); sp = sin( phi( kevt )*deg2rad(1) );
  ct = cos( theta( kevt )*deg2rad(1) ); st = sin( theta( kevt )*deg2rad(1) );
  u  = [ cp*st, sp*st, ct ];
  U  = ones( Na, 1 )*u;
  Di = sqrt( sum( dX.^2, 2 ) - sum( dX.*U, 2 ).^2 );
  Ai = s0( kevt ) - 20/log( 10 )/lambda( kevt )*Di;
  plot( Di, Ai, 'k-' );
  hold on;
  plot( Di, Sa( :, kevt ), 'ks', 'MarkerFace', 'k' );
  hold off; grid on;
  xlabel( 'Lateral distance  [ m ]', 'FontSize', 18 );
  ylabel( 'Amplitude  [ dB ]', 'FontSize', 18 );
  set( gca, 'FontSize', 16 );
  
  Nb   = 100; %% Projection in z = z0 plan
  xmin = -150; xmax = 0;
  dx   = ( xmax - xmin )/Nb;
  xb   = [ xmin:dx:xmax ];
  ymin = -150; ymax = 150;
  dy   = ( ymax - ymin )/Nb;
  yb   = [ ymin:dy:ymax ];
  Xb   = ones( size( yb' ) )*( xb - x0( kevt ) );
  Yb   = ( yb' - y0( kevt ) )*ones( size( xb ) );
  RHO  = sqrt( Xb.^2*( 1 - u( 1 )^2 ) + Yb.^2*( 1 - u( 2 )^2 ) );
  subplot( 1, 2, 2 );
  contour( xb, yb, s0( kevt ) - 20/log( 10 )/lambda( kevt )*RHO, 20 );
  hold on;
  plot( Xa( :, 1 ), Xa( :, 2 ), 'ks', 'MarkerFace', 'k', 'MarkerSize', 12 );
  plot( x0( kevt ), y0( kevt ), 'rp', 'MarkerFace', 'r', 'MarkerSize', 12 );
  plot( X0True( 1 ), X0True( 2 ), 'g^', 'MarkerFace', 'g', 'MarkerSize', 12 );
  hold off; grid on;
  xlabel( 'X  [ m ]', 'FontSize', 18 );
  ylabel( 'Y  [ m ]', 'FontSize', 18 );
  set( gca, 'FontSize', 16 );
  axis equal;
  colorbar;

function myHist( X, nbin, clr, label )

  sX = std( X );
  K = find( label ~= '\' );
  fprintf( [ label( K ), ' = %12.5e\n' ], sX );

  [ N, C ] = hist( X, nbin );
  dC = mean( diff( C ) );
  nrm = sum( N )*dC;
  dN = sqrt( N )/nrm;
  N  = N/nrm;
  h = errorbar( C, N, dN, [ clr, 's' ] );
  set( h( 2 ), 'MarkerFace', clr );
  grid on;
  xlabel( label, 'fontSize', 18 );
  set( gca, 'fontSize', 16 );