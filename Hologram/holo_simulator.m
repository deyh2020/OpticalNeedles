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



% Retuns the ElectricField distribution, as an images for each complex
% component, at the conjugate plane of C_x and C_y, whose are determined by
% the pointres gl1 and gl2 on SLMs
% [C_img,C_img]=f([img],[img]) -->  [a+bi,c+di]=f([1 256],[1 256])
%
function [E_x,E_y] = holo_simulator(slm1,slm2,NoiseAmp,NoisePhase)

N    = size(slm1);

C_SLM1 = zeros(N);
C_SLM2 = zeros(N);

%% Load_SLMresponse;
mapa1 = dlmread('response_SLM1.txt','',2,0);
ph_slm1 = mapa1(:,3);
T_slm1  = mapa1(:,2);

mapa2 = dlmread('response_SLM2.txt','',2,0);
ph_slm2 = mapa2(:,3);
T_slm2  = mapa2(:,2);

%% Exact response: from gray level to complex transmittance
C_SLM1(:) = T_slm1(slm1(:)) .* exp( 1i*ph_slm1(slm1(:)) );
C_SLM2(:) = T_slm2(slm2(:)) .* exp( 1i*ph_slm2(slm2(:)) );


%% Noising the response

% white noise mask to be added to the response 
%   MASK = MAX + (MAX - MIN).*rand(N)
mask1_A = NoiseAmp(1)   + (NoiseAmp(2) - NoiseAmp(1)  ).*rand(N);
mask1_P = NoisePhase(1) + (NoisePhase(2)-NoisePhase(1)).*rand(N);
mask2_A = NoiseAmp(1)   + (NoiseAmp(2) - NoiseAmp(1)  ).*rand(N);
mask2_P = NoisePhase(1) + (NoisePhase(2)-NoisePhase(1)).*rand(N);


C_SLM1 = C_SLM1.*exp(1i*mask1_P)+mask1_A;
C_SLM2 = C_SLM2.*exp(1i*mask2_P)+mask2_A;

% figure
% imagesc(abs2(C_SLM1));  title |E^{SLM}_x|^2
% figure
% imagesc(angle(C_SLM1)); title \phi^{SLM}_x 
% figure
% imagesc(abs2(C_SLM2));  title |E^{SLM}_y|^2
% figure
% imagesc(angle(C_SLM2)); title \phi^{SLM}_y


%% Arizon procedure

FC_x = fftWELL( C_SLM1 );
FC_y = fftWELL( C_SLM2 );

spatial_filter = zeros(N);
spatial_filter( ceil(3*N(1)/8):ceil(5*N(1)/8)-1 , ...
	            ceil(3*N(2)/8):ceil(5*N(2)/8)-1 )  = 1;

FC_x = FC_x.*spatial_filter;
FC_y = FC_y.*spatial_filter;

E_x = ifftWELL( FC_x );
E_y = ifftWELL( FC_y );

% figure
% imagesc(abs2(E_x));   title |E_x|^2
% figure
% imagesc(abs2(E_y));   title |E_y|^2
% figure
% imagesc(angle(E_x)); title \phi_x
% figure
% imagesc(angle(E_y)); title \phi_y
% figure
% imagesc(abs2(E_x)+abs2(E_y)); title |E|^2


