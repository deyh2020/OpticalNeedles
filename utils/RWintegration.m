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


% Integrates the Richards-Wolf spectrum of plane wave
% where E0x and E0y are the X nad Y components of the spectrum.
% E0x and E0y must be 2D-matrix where the first dim is the THETA coordinate
% and the second dim is the PHI coordinate.

clear variables;
E0x=ones(512,513);E0y=ones(512,513);theta0=60*pi/180;z=1;r_max=10;RES=128;  


%function[EFx,EFy,EFz,xd,yd,u,v]=RWintegration(E0x,E0y,theta0,z,r_max,RES)

% size of matrix
N    = size(E0x);
oneT = ones(1,N(1));
oneP = ones(1,N(2));

% Pupils coordinates
theta = linspace(0,theta0,N(1));
phi   = linspace(0, 2*pi ,N(2));

% Focusing zone coordinates
if max(size(RES))==1,RES(2)=RES(1);end
r   = linspace(0,r_max,RES(1));
PHI = linspace(0, 2*pi,RES(2));


% azimuthal and radial vectors
%  e1 = [ -sin(phi) , cos(phi) , 0 ]
e1x =-oneT.'*sin(phi);
e1y = oneT.'*cos(phi);
e1z = zeros(N);
%  e2 = [ cos(theta)cos(phi) , cos(theta)sin(phi) , sin(theta) ]
e2x = cos(theta).'*cos(phi);
e2y = cos(theta).'*sin(phi);
e2z = sin(theta).'*oneP;

% azimuthal and radial components
f1 = E0x.*e1x + E0y.*e1y ; % f1 = E0 · e1
f2 = E0x.*e2x + E0y.*e2y ; % f2 = E0 · e2

for k=1:length(z)
zk=z(k);
    for i=1:RES(1)
    ri=r(i);
        for j=1:RES(2)
        PHIj=PHI(j);
            
            ex = @(varP,varT) (
        
            Ex(i,j,k) = integral2( 
        
        
        
        end
    end
end








