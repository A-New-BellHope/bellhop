% Parabolic Bottom test
% mbp

global units jkpsflag
units = 'km';

% *************************************************
% linear boundary interpolation

make_bdry( 'L' )

%%
% the rays:
bellhop ParaBot
figure; plotray ParaBot
axis( [ 0 20 -5000 5000 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGeom.ati' );

bellhop ParaBotTLGeom

plotshd( 'ParaBotTLGeom.shd', 1, 1, 1 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

%%
% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGB.ati' );

bellhop ParaBotTLGB

plotshd( 'ParaBotTLGB.shd', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

%print -dpng ParabotLTL
%%

% *************************************************
% curvilinear boundary interpolation

make_bdry( 'C' )

%%
% the rays:
bellhop ParaBot
figure; plotray ParaBot
axis( [ 0 20 -5000 5000 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGeom.ati' );
bellhop ParaBotTLGeom
plotshd( 'ParaBotTLGeom.shd', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGB.ati' );
bellhop ParaBotTLGB
plotshd( 'ParaBotTLGB.shd', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot
