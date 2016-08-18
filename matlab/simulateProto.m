clear all;
close all;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
%%                             Antenna Settings
%%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xa = [ -2.510000e+001   1.305000e+002   2.663200e+003;
       -5.210000e+001   8.360000e+001   2.664400e+003;
       -2.340000e+001   5.100000e+000   2.671200e+003;
       -1.144000e+002  -2.380000e+001   2.669600e+003;
       -4.700000e+001  -7.370000e+001   2.670100e+003;
       -2.160000e+001  -1.104000e+002   2.670800e+003  ];
Na = size( Xa, 1 );
FS = 1/5e-9;
C0 = 3e8;

Nevt = 1000;
dT0  = 0.5;  %% Error in samples
dA0  = 0.05; %% Relative error on amplitudes / linear scale

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
%%                              Point Source
%%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Settings
%%===
Xs = mean( Xa, 1 );
Ts = 0;
CR = 0.9;

%% Arrival times: point source like
%%===
Da = sqrt( sum( ( Xa - ones( Na, 1 )*Xs ).^2, 2 ) );
Ta = Da/( C0*CR ) + Ts;

%% Amplitude: 1/r attenuation
%%===
S0 = 20*log( 1024 );
Ai = S0 - 20*log10( Da );

%% Simulate data
%%===
dT   = dT0/FS*randn( Na, Nevt );
dA   = dA0*randn( Na, Nevt );
fid = fopen( 'pointsource-data.txt', 'w+' );
for k = 1:Nevt
  Tm = Ta + dT( :, k );
  Tm = Tm - min( Tm );
  Am = Ai.*( 1 + dA( :, k ) );
  fprintf( fid, ' %4d', Na );
  for ka = 1:Na
    fprintf( fid, ' %12.6e', [ Xa( ka, : ), Tm( ka )*C0, Am( ka ) ] );      
  end
  fprintf( fid, '\n' );
end
fclose( fid );

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
%%                    Plan Wave & Exponentiel Amplitude Model
%%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Settings
%%===
theta  = 20.0*deg2rad(1);
phi    = 60*deg2rad(1);
u      = [ cos( phi )*sin( theta ), sin( phi )*sin( theta ), cos( theta ) ];
CR     = 0.9;

X0     = mean( Xa, 1 );
S0     = 20*log10( 1024 );
lambda = 100;
a0     = 20/log( 10 )/lambda;

%% Arrival times: plan wave like
%%===
dX = [ Xa - ones( Na, 1 )*X0 ];
U  = ones( Na, 1 )*u;
Ta = sum( dX.*U, 2 )/( C0*CR );

%% Amplitudes: exponentiel model
%%===
Di = sqrt( sum( dX.^2, 2 ) - sum( dX.*U, 2 ).^2 );
Ai =  S0 - a0*Di;

%% Simulate data
%%===
dT   = dT0/FS*randn( Na, Nevt );
dA   = dA0*randn( Na, Nevt );
fid = fopen( 'planwave-data.txt', 'w+' );
for k = 1:Nevt
  Tm = Ta + dT( :, k );
  Tm = Tm - min( Tm );
  Am = Ai.*( 1 + dA( :, k ) );
  fprintf( fid, ' %4d', Na );
  for ka = 1:Na
    fprintf( fid, ' %12.6e', [ Xa( ka, : ), Tm( ka )*C0, Am( ka ) ] );      
  end
  fprintf( fid, '\n' );
end
fclose( fid );