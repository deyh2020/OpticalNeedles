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



% Compute the high focused field in the focal zone by means of
% Eq. (3.47) of Novotny's book.
% Where the integration is done via Fast Fourier Transform function.
%
% Einc1, Einc2 are the complex amplitudes at the Entrance Pupil (paraxial)
% with a certain base specified in the last argument: 
%   - 'Cartesian': (Einc1=Ex and Einc2=Ey)
%   - 'AzimuRadial': base (Einc1=E0·e1 and Einc2=E0·e'2)
% 
% The fisical width of the inpust must be [-N/2/L (N/2/L-1/L] 
% where N=size(Einc1).
%    ***    It can be useful to compute the Einc using the    ***
%    ***   coordinates returned by coord2RWinFFT() function.  ***
%
% NA is the Numerical Aperture of an isoplanatic microscope objective.
% L is the full-width of the focal window.
% Z is the distance from the focal plane.
%    *** Z can be an array to compute a 3D-focal zone.  ***
%
function[EFx,EFy,EFz]=RWinFFT(Einc1,Einc2,NA,N,L,Z,base)

if N < 2*NA^2*max(abs(Z))/sqrt(1-NA^2),warning('N is in undersampling conditions!');end

%% Coordinates on the Pupil

[theta,phi,mask] = coord2RWinFFT(N,L,NA);

% spherical coordinates on the exit pupil
sinfi = sin(phi)   .*mask;
cosfi = cos(phi)   .*mask;
sinte = sin(theta) .*mask;
coste = cos(theta) .*mask;


%% Vectors and parameters of Richards-Wolf integral

% appodization of isoplanatic an objective microscope
P = sqrt(coste); 

% e1 vector (Azimuthal)
e1x = - sinfi.*mask;  
e1y =   cosfi.*mask;
e1z =   0;

% e2 vector (Radial on the Gaussian sphere)
e2x =   coste.*cosfi.*mask;
e2y =   coste.*sinfi.*mask;
e2z = - sinte.*mask;

% ep vector (Radial on the Paraxial beam)
epx =   cosfi.*mask;
epy =   sinfi.*mask;
epz =   0;



%% Angular spectrum of plane waves derivation

switch base
    case 'Cartesian'
        % azimuthal and radial components
        f1 = Einc1.*e1x + Einc2.*e1y ; % f1 = E0 · e1
        f2 = Einc1.*epx + Einc2.*epy ; % f2 = E0 · e'2
    case 'AzimuRadial'
        f1 = Einc1;
        f2 = Einc2;
    otherwise
        error('Choose a valid base ("Cartesian" or "AzimuRadial")')
end

EFx = zeros(N,N,length(Z));
EFy = zeros(N,N,length(Z));
EFz = zeros(N,N,length(Z));
wb = waitbar( 0, 'Richards-Wolf integration via FFT. Z-loop:' );
for iz=1:length(Z)
z = Z(iz); % current z
    % Angular spectrum of plane waves
    Vx = P.*( f1.*e1x + f2.*e2x ) .*exp( 1i*2*pi*coste.*z );
    Vy = P.*( f1.*e1y + f2.*e2y ) .*exp( 1i*2*pi*coste.*z );
    Vz = P.*( f1.*e1z + f2.*e2z ) .*exp( 1i*2*pi*coste.*z );

    RWfact = 1;% -1i/N/N; % /lambda/f (to calculate the real energy)
    
    EFx(:,:,iz) = fftshift(fft2(ifftshift( Vx./(coste+eps).*mask )))*RWfact; 
    EFy(:,:,iz) = fftshift(fft2(ifftshift( Vy./(coste+eps).*mask )))*RWfact;
    EFz(:,:,iz) = fftshift(fft2(ifftshift( Vz./(coste+eps).*mask )))*RWfact;
    
    if length(Z)>1
        waitbar( iz/length(Z) , wb );
    end
end
delete(wb);


