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



% Computes the propagated field following Eq. (3.47) of Novotny's book.
% The integration is via trapz() function
% [E0x,E0y] are the complex distributions on the entrance pupil (paraxial)
% [xx,yy,zz] are the array in the focal zone where field is calculated
% NA is the Numerical Aperture of an isoplanatic microscope objective

%---F-O-R---T-E-S-T-I-N-G-----
% E0x = ones(200,200);       %
% E0y = ones(200,200);       %
% NA  = 0.85;                %
% xx  = linspace(-4,4,50);   %
% yy  = linspace(-4,4,50);   %
% zz  = 0;                   %
%-----------------------------

function[EFx,EFy,EFz]=RWinTrapz(E0x,E0y,NA,xx,yy,zz)

N  = size(E0x,1);

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

P     = sqrt(coste); % APPODIZATION FACTOR !!! P = sqrt(coste) !!!
 
% vector e1
e1x = - sinfi.*mask;           
e1y =   cosfi.*mask;
e1z =   0;

% vector e2
e2x =   coste.*cosfi.*mask;      
e2y =   coste.*sinfi.*mask;
e2z = - sinte.*mask;

% vector e'2
epx =   cosfi.*mask;      
epy =   sinfi.*mask;
epz =   0;

% azimuthal and radial components
f1 = E0x.*e1x + E0y.*e1y ; % f1 = E0 · e1
f2 = E0x.*epx + E0y.*epy ; % f2 = E0 · e'2

% initializing the Focused beam
EFx = zeros(length(yy),length(xx),length(zz));
EFy = zeros(length(yy),length(xx),length(zz));
EFz = zeros(length(yy),length(xx),length(zz));
for k = 1:length(zz)
z=zz(k);

    % appodized Plane Wave Angular Spectrum in propagation
    %        sqrt(theta)·E0·exp{-ikzcos(theta)}
    Vx = P.*( f1.*e1x + f2.*e2x ).*exp( -1i*2*pi*coste.*z );
    Vy = P.*( f1.*e1y + f2.*e2y ).*exp( -1i*2*pi*coste.*z );
    Vz = P.*( f1.*e1z + f2.*e2z ).*exp( -1i*2*pi*coste.*z );
    
    for i = 1:length(xx)
    xi=xx(i);
        for j=1:length(yy)
        yi=yy(j);

            rr  = sqrt(xi.^2+yi.^2);
            PHI = atan2(yi,xi);

            Fact = exp( 1i*2*pi * rr*sinte .* cos(PHI-fi) );
                
               % doble integral via TRAPZ(y, TRAPZ(x,f(x,y)) )
            EFx(j,i,k) = trapz(y(:,1) ,  trapz( x(1,:) , Vx.*Fact ,2)  );                          %calcula FFT2
            EFy(j,i,k) = trapz(y(:,1) ,  trapz( x(1,:) , Vy.*Fact ,2)  );
            EFz(j,i,k) = trapz(y(:,1) ,  trapz( x(1,:) , Vz.*Fact ,2)  ); 

        end
        if N>100
        multiWaitbar( 'Richards-Wolf integration via Trapz. X-loop:', i/length(xx) );
        end
    end
    if length(zz)>1
    multiWaitbar( 'Richards-Wolf integration via Trapz. Z-loop:', k/length(zz) );
    end
end




