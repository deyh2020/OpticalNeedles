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
% Ho fa utilitzant transformadea de Fourier.
% La funció d'entrada és (Eincx,Eincy,Eincz)". 
% Utilitza la funció fhermite_d(D,N,Lm,n,m)per generar un hermite d'ordre n,m. 
% La funció fhermite(D,N,Lm,n,m)genera un hermite d'ordre n,m i fa la transformada per tenir-lo a l'espai de Fourier. 


function[EFx,EFy,EFz]=RWinFFT(f1,f2,z,NA,L)

N=size(f1,1);



[x1,y1] = meshgrid( linspace(-1,1,N) , linspace(-1,1,N) );  %genera matriu


ro2 = x1.*x1+y1.*y1;            %radi a la pupil·la
ro  = sqrt(ro2);
mask= (ro<=NA).*1;
ro  = mask.*ro;

fi    = atan2(y1,x1);
sinfi = sin(fi);
cosfi = cos(fi);
sinte = ro;
coste = sqrt(1-sinte.^2);

P     = sqrt(coste);


%%
 
e1x =-sinfi.*mask;                                          % vector e1
e1y = cosfi.*mask;
e1z = 0;

e2x = coste.*cosfi.*mask;                                  % vector e2
e2y = coste.*sinfi.*mask;
e2z =-sinte.*mask;

% azimuthal and radial components
%f1 = E0x.*e1x + E0y.*e1y ; % f1 = E0 · e1
%f2 = E0x.*e2x + E0y.*e2y ; % f2 = E0 · e2

% plane wave angular spectrum E0
M=32;Mm=floor(M/2)*N+1;MM=(floor(M/2)+1)*N;
Vx = zeros(M*N,M*N) ; Vy = zeros(M*N,M*N) ; Vz = zeros(M*N,M*N);
Vx(Mm:MM,Mm:MM) = P.*( f1.*e1x + f2.*e2x ).*exp( 1i*2*pi*coste.*z );
Vy(Mm:MM,Mm:MM) = P.*( f1.*e1y + f2.*e2y ).*exp( 1i*2*pi*coste.*z );
Vz(Mm:MM,Mm:MM) = P.*( f1.*e1z + f2.*e2z ).*exp( 1i*2*pi*coste.*z );

EFx=fftshift(fft2(ifftshift(Vx)));%/N/M/N/M;                           %calcula FFT2
EFy=fftshift(fft2(ifftshift(Vy)));%/N/M/N/M; 
EFz=fftshift(fft2(ifftshift(Vz)));%/N/M/N/M; 

% TVx=abs2(Vx);                         %calcula el modul quadrat
% TVy=abs2(Vy);
% TVz=abs2(Vz);
% 
% IVx=sum(sum(TVx));                      %calcula l'energia
% IVy=sum(sum(TVy));
% IVz=sum(sum(TVz));
% IV=IVx+IVy+IVz;
% 
% TEFx=abs2(EFx);                         %calcula el modul quadrat
% TEFy=abs2(EFy);
% TEFz=abs2(EFz);
% TEFtrans=TEFx+TEFy;
% TEFTOT=TEFx+TEFy+TEFz;
% 
% IEFx=sum(sum(TEFx));                      %calcula l'energia
% IEFy=sum(sum(TEFy));
% IEFz=sum(sum(TEFz));
% IEFtrans=sum(sum(TEFtrans));
% IEFTOT=sum(sum(TEFTOT));
% IEFt=IEFx+IEFy;
% IEF=IEFx+IEFy+IEFz;



