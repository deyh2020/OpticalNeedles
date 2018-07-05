%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Applied Physics and Optics (UB) and Optics and Laser Beams (UCM)
% 
%                 https://github.com/dmaluenda/OpticalNeedles
%
%                 David Maluenda Niubo - dmaluendn@gmail.com            
%
%      CC: by, NC, SA                                         2012-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Returns the theta and phi coordinates and an Entrance Pupil mask
% of size N and with a total physical size of L (in lambdas) and 
% Numerical Aperture NA.
%
% ** This Function is useful to be employed with the RWinFFT() function **
%
function [theta,phi,mask] = coord2RWinFFT(N,L,NA)

% x1 and y1 are the coordinates on the Entrance Pupil in lambda^{-1} units
[x1,y1] = meshgrid( linspace(-N/2,N/2,N)/L , linspace(N/2,-N/2,N)/L );

% radial coordinate on the entrance pupil
rho2 = x1.*x1 + y1.*y1;
rho  = sqrt(rho2);
mask = (rho<=NA).*1;

theta = asin(rho);
phi   = atan2(y1,x1);