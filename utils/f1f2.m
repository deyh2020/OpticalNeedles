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



% Reurns the Azimutal f1 and Radial f2 transversal components of a
% Cartessian field Ex and Ey with a certain Numerical Aperture NA

%---F-O-R---T-E-S-T-I-N-G-----
% E0x = ones(200,200);       %
% E0y = ones(200,200);       %
% NA  = 0.85;                %
% zz  = 0;                   %
% xx  = linspace(-4,4,50);   %
% yy  = linspace(-4,4,50);   %
%-----------------------------

function[f1,f2]=f1f2(Ex,Ey,NA)

N  = size(Ex,1);

% Exit Pupil variables
[x,y] = meshgrid( linspace(-NA,NA,N) , linspace(NA,-NA,N) );

ro2 = x.*x+y.*y;   
ro  = sqrt(ro2);
mask= (ro<=NA).*1;
ro  = mask.*ro;

fi    = atan2(y,x);
sinfi = sin(fi);
cosfi = cos(fi);
sinte = ro;
coste = sqrt(1-sinte.^2);

P     = sqrt(coste);
 
% vector e1
e1x = -sinfi.*mask;           
e1y =  cosfi.*mask;
e1z =  0;

% vector e2
e2x =  coste.*cosfi.*mask;      
e2y =  coste.*sinfi.*mask;
e2z = -sinte.*mask;

% azimuthal and radial components
f1 = Ex.*e1x + Ey.*e1y ; % f1 = E0 · e1
f2 = Ex.*e2x + Ey.*e2y ; % f2 = E0 · e2





