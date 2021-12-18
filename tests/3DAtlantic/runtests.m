% 3D Normal mode tests using N. Atlantic scenario
% mbp

global units
units = 'km';

% plot the triangulation

figure;
plottri( 'lant' )
axis( [ 100 600 200 600 ] )
axis( [ 0 650 0 650 ] )
%axis( [ 50 600 150 600 ] )

% print -depsc2 fig5_19b.eps

%%

krakenall

%%

field3d lant

% field3dM is broken
% Looks like it has not been updated to fully handle the freqVec now used
% for broadband

%field3dM lant

figure
plotshdpol( 'lant.shd', 333, 315, 400 )
axis( [ 0 650 0 650 ] )
caxisrev( [ 80 100 ] )

% Gaussian beam version
% eval( [ '! "' runfield3d '" lant_gbt > field3d.prt' ] );

% figure;
% plotshdpol 'lant_gbt.shd'
% title( 'F = 50 Hz, SD = 400 m, RD = 400 m' )
% axis( [ 0 650 0 650 ] )
% caxisrev( [ 80 100 ] )

% print -depsc2 fig5_20.eps

%%

% test of a single bearing calculation
% need to change .env files to include multiple depths first

field3d lant_90degrees

figure
plotshd 'lant_90degrees.shd'
title( 'F = 50 Hz, SD = 400 m, theta = 90 degrees' )
caxisrev( [ 70 100 ] )


% test of a volume calculation

%!field3d < lantvol.flp
%copyfile 'SHDFIL' 'lantvol.shd'
%plotshdpol 'lantvol.shd'

%axis( [ 100 600 200 600 ] )
%view( 0, -90 )

% Gaussian beam run (takes several hours)
%!field3d.exe < lantgbt.flp
%movefile 'SHDFIL' 'lantgbt.shd'
%plotshdpol 'lantgbt.shd'
%title( 'F = 50 Hz, SD = 400 m, RD = 400 m' )
