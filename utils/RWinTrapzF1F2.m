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


% Calcula el feix propagant en el cas convergent segons la formula (3.47) del Novotny.
% Integra via la funció Trapz per a 
% [E0x,E0y] és la distribució complexa de les components X i Y en la pupila d'entrada.
% [xx,yy,zz] són els arrays de coordenades en la zona focal

%---F-O-R---T-E-S-T-I-N-G-----
% E0x = ones(200,200);       %
% E0y = ones(200,200);       %
% NA  = 0.85;                %
% zz  = 0;                   %
% xx  = linspace(-4,4,50);   %
% yy  = linspace(-4,4,50);   %
%-----------------------------

function[EFx,EFy,EFz]=RWinTrapzF1F2(E01,E02,NA,xx,yy,zz)

N  = size(E01,1);

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
e1x = -sinfi.*mask;           
e1y =  cosfi.*mask;
e1z =  0;

% vector e2
e2x =  coste.*cosfi.*mask;      
e2y =  coste.*sinfi.*mask;
e2z = -sinte.*mask;

% azimuthal and radial components
f1 = E01;
f2 = E02;

% initializing the Focused beam
EFx = zeros(length(yy),length(xx),length(zz));
EFy = zeros(length(yy),length(xx),length(zz));
EFz = zeros(length(yy),length(xx),length(zz));
if length(zz)>1
multiWaitbar( 'Richards-Wolf integration via Trapz. Z-loop:', 0, 'Color', 'g', 'CanCancel', 'on' );
end
for k = 1:length(zz)
z=zz(k);

    % appodized Plane Wave Angular Spectrum in propagation
    %        sqrt(theta)·E0·exp{-ikzcos(theta)}
    Vx = P.*( f1.*e1x + f2.*e2x ).*exp( -1i*2*pi*coste.*z );
    Vy = P.*( f1.*e1y + f2.*e2y ).*exp( -1i*2*pi*coste.*z );
    Vz = P.*( f1.*e1z + f2.*e2z ).*exp( -1i*2*pi*coste.*z );

    multiWaitbar( 'Richards-Wolf integration via Trapz. X-loop:', 0, 'Color', 'y', 'CanCancel', 'on' );
    for i = 1:length(xx)
    xi=xx(i);
        for j=1:length(yy)
        yi=yy(j);

            rr  = sqrt(xi.^2+yi.^2);
            PHI = atan2(yi,xi);

            Fact = exp( 1i*2*pi * rr*sinte .* cos(PHI-fi) );
                
               % doble integral via TRAPZ(y, TRAPZ(x,f(x,y)) )
            EFx(j,i,k) = trapz(y(:,1) ,  trapz( x(1,:) , Vx.*Fact,2)  );                          %calcula FFT2
            EFy(j,i,k) = trapz(y(:,1) ,  trapz( x(1,:) , Vy.*Fact,2)  );
            EFz(j,i,k) = trapz(y(:,1) ,  trapz( x(1,:) , Vz.*Fact,2)  ); 

        end
        multiWaitbar( 'Richards-Wolf integration via Trapz. X-loop:', i/length(xx) );
    end
    if length(zz)>1
    multiWaitbar( 'Richards-Wolf integration via Trapz. Z-loop:', k/length(zz) );
    end
end
multiWaitbar( 'CloseAll' );


% 
% figure;imagesc(xx,yy,abs(EFx))
% figure;imagesc(xx,yy,abs(EFy))
% figure;imagesc(xx,yy,abs(EFz))
% figure;imagesc(xx,yy,abs2(EFx)+abs2(EFy)+abs2(EFz))
% 
% TVx=abs2(Vx);                         %calcula el modul quadrat
% TVy=abs2(Vy);
% TVz=abs2(Vz);
% 
% IVx=sum(sum(TVx))                      %calcula l'energia
% IVy=sum(sum(TVy))
% IVz=sum(sum(TVz))
% IV=IVx+IVy+IVz


