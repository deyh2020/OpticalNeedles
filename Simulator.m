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


clear variables; 
close all;
tic; 


%% Parameters and Variables
verbose = 0;  % figures displaying mode -> 0: normal ; 1: detailed ; 2: extradetailed

% design of the AngularSpectrum
NA        = 0.65; 
Npix      = 180; % Diameter of Entrance pupil in SLM's pixels
ill_FLAG  = 'Scalar'; % Choose either 'Vectorial' or 'Scalar' for the illumination
doHOLO    = 1;  % 0: loads the holograms from 'Hologram<SLMnumber>.bmp'
                % 1: computes the hologram following the next parameters
BEAM_flag = 1;  % 0: h(a)=1 ; 1: h(a)=sinc(a) ; 2: h(a)=SUM[sin(a)] ; 3: h(a)=N·exp[-(N alpha)^2]
Nbeam     = 8;  % N parameter Eq.(11) -> h(alpha) = N sinc(2 pi N ...) 
mbeam     = 8;  % m parameter Eq.(12) -> alpha(m) = m(1-alpha_0)/(2N)+alpha_0
doSphAbCorr = 0; % Correct the Spherical Aberration? 0->NO ; any->YES 


% simulation of the optical system
NOISE_flag = 1;  % CHOOSE: 1=noisy simulation ; 0=unnoisy simulation 
% noise from errors in the SLM characterization 
wNoise_Amp = [-0.2 0.2];  % white noise factor for SLM's amplitude characterization (exp. +/-2%)
wNoise_Ph  = [-30 30]*pi/180; % white noise factor for SLM's phase characterization (exp. +/-10º)
% noise from other error sources like polarizers or almost-planar mirrors
pNoise_Amp = [0.1 0.01];  % Perlin noise as amplitude backgroung from polarizers (exp. 5% -both positive-)
pNoise_Ph  = [-30 30]*pi/180; % Perlin noise as for phase background from mirrors distortion (exp. +/-36º=lambda/10)
% aberrations and optical alignment
As       = -10;     % spherical aberration's coefficient in lambdas (exp. As=-10)
disAlign = [0 0]; % descentering the hologram from the optical axis [Dx Dy] 


% Richards-Wolf integration
thetaM = asin(NA);  % max of theta RW-integration, i.e. asin(NA)
f0     = 1.5;       % Filling Factor of Gaussian Beam [f0->Infty => PlaneWave]
N      = 769;       % Pixel-resolution of the FFT procedure
Dz     = [-50 50];  % [Zmin Zmax] in lambdas 
Nz     = 101;       % Pixel-resolution of the Z-coordinate
Lxy    = 12;        % XY half-size ROI in Focal Zone
X_0    = 0;         % Certain value for each coordinate
Y_0    = 0;         %    when a constant is recquired such as for
Z_0    = 0;         %    plotting the information (cuts and profiles)

if ~doSphAbCorr
    Dz=Dz+As;
end



%% Variables and Coordinates

% preparing the Richards-Wolf integration via FFT
L       = 130;%round(Npix*3/4/2)*2; % full-size of FFTwindow (in lambdas) / Resolution of Entrance Pupil
xy      = linspace(-L/2,L/2,N);            % array of x positions (or y);
if Nz > 1, Zs = linspace(Dz(1),Dz(2),Nz);  % array of z positions
else Zs = Z_0; end % if Nz is 1 the simulations is a cut on Z_0

% getting the coordinates on the Pupil
[theta,phi,mask] = coord2RWinFFT(N,L,NA);

% spherical coordinates on the exit pupil
sinfi = sin(phi)    .*mask;
cosfi = cos(phi)    .*mask;
sinte = sin(theta)  .*mask;
coste = cos(theta)  .*mask;
rho   = asin(theta) .*mask;

% Regions Of Interest [just for plotting]
% ROI in focal zone
[~,pm] = min(abs( xy+Lxy )); % index of -Lxy value
[~,pM] = min(abs( xy-Lxy )); % index of +Lxy value
XY     = xy(pm:pM);          % XY-coordinate array (focal zone)
% Nxy    = pM-pm+1;            % pixel resolution of XY coordinate array
[~,pX] = min(abs(XY-X_0));   % index of X_0
[~,pY] = min(abs(XY-Y_0));   % index of Y_0
[~,pZ] = min(abs(Zs-Z_0));   % index of Z_0
% ROI in paraxial zone
SLM_ROI = ceil(N/2-Npix/2) : ceil(N/2+Npix/2);



%% Incident beam

% Polarization of the illumination
if strcmp(ill_FLAG,'Vectorial')
    ill_X = cosfi;
    ill_Y = sinfi;
elseif strcmp(ill_FLAG,'Scalar')
    ill_X = ones(size(cosfi));
    ill_Y = zeros(size(sinfi));
else
    error('Choose a correct illumination condition (Vectorial or Scalar)');
end

% Aberrated Gaussian input beam
g = exp(-(sinte./sin(thetaM)/f0).^2).*exp(1i*2*pi*As*sinte.^4).*mask; 


if doHOLO 
% Compute the hologram 
    alpha  = cos(theta).*mask;
    alpha0 = cos(thetaM);
    alphaB = mbeam/2/Nbeam*(1-alpha0) + alpha0 ; % alpha(m) in the paper

    if Nbeam==0
        fn = 1;
    else
        % different modulations types
        switch BEAM_flag
            case 1  % f(alpha) = N·sinc(...) -> Paper:Eq.(11) 
                fn = Nbeam * sinc( 2*Nbeam * (alpha-alphaB)./(1-alpha0) ) ;

            case 2  % f(alpha) = SUM{ sin(...)·sign(*) } -> first approx on March'16
                Ns = 1:Nbeam;
                fn = zeros(N);
                for i=1:length(Ns)
                    Ni=Ns(i);
                    fi = sin( pi*Ni * (alpha-alpha0)./(1-alpha0) ) .* ...
                         sign(real((1i+.01).^(Ni+1))); %  -1   1  -1   1  -1   1  -1

                    fn = fn+fi;
                end

            case 3  % f(a)=N·exp[-(N alpha)^2]
                fn = Nbeam/sqrt(pi)*exp(-Nbeam.^2 * (alpha-alphaB).^2 ./ (1-alpha0).^2 );

            otherwise  % f(a)=1   -case 0: DEFAULT-
                fn = 1;
        end
    end

    holo = fn.*mask.*sinte;    % h(theta) = f(theta) * sin(theta)
    holo(isnan(holo)) = 0;
    holo = holo/max(holo(:));

    if doSphAbCorr
        cor_sph = exp(-1i*2*pi*As*sinte.^4);
    else
        cor_sph = 1;
    end
    
    % Designed distribution
    EDx = holo.*cor_sph.*mask.*ill_X;
    EDy = holo.*cor_sph.*mask.*ill_Y;

    %% Plot design
    if verbose > 0
        % the Azimuthal and Radial components from the cartesians
        [f1_D,f2_D] = f1f2( EDx , EDy , NA );

        % Designed beam
        f1_D = f1_D.*g.*mask;
        f2_D = f2_D.*g.*mask;

        FontSize = 18;
        figure;
        
        ax=subplot(2,2,1); imagesc(abs2(EDx(SLM_ROI,SLM_ROI)),[0 1]);
        cmap('hot'); title('$|E_D \cdot e_x|^2$','Interpreter','Latex') ; axis square
        ax.FontSize = FontSize;
        ax.LineWidth = 2;
        ax.XTick = 0:60:180;
        ax.YTick = 0:60:180;
        axis([0 180 0 180])
        
        ax=subplot(2,2,2); imagesc(abs2(EDy(SLM_ROI,SLM_ROI)),[0 1]);
        cmap('hot'); title('$|E_D \cdot e_y|^2$','Interpreter','Latex')  ; axis  square
        ax.FontSize = FontSize;
        ax.LineWidth = 2;
        ax.XTick = 0:60:180;
        ax.YTick = 0:60:180;
        axis([0 180 0 180])
        
        ax=subplot(2,2,3); imagesc(abs2(f1_D(SLM_ROI,SLM_ROI)),[0 1]);
        cmap('hot'); title('$|E_D \cdot e_1|^2$','Interpreter','Latex')  ; axis  square
        ax.FontSize = FontSize;
        ax.LineWidth = 2;
        ax.XTick = 0:60:180;
        ax.YTick = 0:60:180;
        axis([0 180 0 180])
        
        ax=subplot(2,2,4); imagesc(abs2(f2_D(SLM_ROI,SLM_ROI)),[0 1]);
        cmap('hot'); title('$|E_D \cdot e_2|^2$','Interpreter','Latex')  ; axis  square
        ax.FontSize = FontSize;
        ax.LineWidth = 2;
        ax.XTick = 0:60:180;
        ax.YTick = 0:60:180;
        axis([0 180 0 180])
        
        %% in very detailed mode, we show the phases and
        %     compute the exact solution (without the experimental parameters)
        if verbose > 1

            figure; 
            ax=subplot(3,2,1);
            imagesc(angle(EDx(SLM_ROI,SLM_ROI)),[-pi pi]);
            title('$\phi_x$','Interpreter','Latex') ;axis square
            cmap('phase');
            ax.FontSize = FontSize;
            ax.LineWidth = 2;
            ax.XTick = 0:60:180;
            ax.YTick = 0:60:180;
            axis([0 180 0 180])
        
            ax=subplot(3,2,2);
            imagesc(angle(EDy(SLM_ROI,SLM_ROI)),[-pi pi]);
            title('$\phi_y$','Interpreter','Latex');axis square
            cmap('phase');
            ax.FontSize = FontSize;
            ax.LineWidth = 2;
            ax.XTick = 0:60:180;
            ax.YTick = 0:60:180;
            axis([0 180 0 180])
            
            ax=subplot(3,2,[3 4 5 6]);
            imagesc(mod(angle(EDy(SLM_ROI,SLM_ROI))-angle(EDx(SLM_ROI,SLM_ROI)),2*pi),[-pi pi]);
            title('$\phi_y - \phi_x$','Interpreter','Latex') ;
            cmap('phase');
            axis square
            ax.FontSize = FontSize;
            ax.LineWidth = 2;
            ax.XTick = 0:60:180;
            ax.YTick = 0:60:180;
            axis([0 180 0 180])
        
            % Richards-Wolf integration
            [EFx,EFy,EFz] = RWinFFT(f1_D,f2_D,NA,N,L,Zs,'AzimuRadial');

            IF_D  = squeeze( abs2(EFx(pm:pM,pm:pM,:)) + ...
                             abs2(EFy(pm:pM,pm:pM,:)) + ...
                             abs2(EFz(pm:pM,pm:pM,:)) );

            IF_Dt = squeeze( abs2(EFx(pm:pM,pm:pM,:)) + ...
                             abs2(EFy(pm:pM,pm:pM,:)) );

            M_D   = max(sqrt(IF_D(:)));
            IF_D  = IF_D/M_D/M_D;
            IF_Dt = IF_Dt/M_D/M_D;
            EFx_D = EFx(pm:pM,pm:pM,:) / M_D;
            EFy_D = EFy(pm:pM,pm:pM,:) / M_D;
            EFz_D = EFz(pm:pM,pm:pM,:) / M_D;

            clear EFx EFy EFz

            figure;

            hs(1) = subplot(3,1,1);
            plot(XY,squeeze(IF_Dt(pX,:,pZ)),'-','Linewidth',3);hold on
            plot(XY,squeeze(IF_D(pX,:,pZ)),'--','Linewidth',2);
            xlabel('$x$ in $\lambda$','FontSize',20,'Interpreter','latex');
            ylabel('$I(x,y=0,z=0)$','FontSize',20,'Interpreter','latex')
            axis([min(XY) max(XY) 0 1.05]) ; legend 'Transversal' 'Total'
            title(['EXACT solution for: $N = ' num2str(Nbeam) ...
                                   '$ ; $m = ' num2str(mbeam) '$'], ...
                    'Interpreter','latex' , 'FontSize',24 )

            hs(2) = subplot(3,1,2);
            surf(Zs,XY,squeeze(IF_D(:,pY,:)));cmap('red');
            shading interp;axis equal;view(2);
            xlabel('$z$ in $\lambda$','FontSize',20,'Interpreter','latex');
            ylabel('$x$ in $\lambda$','FontSize',20,'Interpreter','latex')
            hs(2).LineWidth = 2;
            
            hs(3) = subplot(3,1,3);
            plot(Zs,squeeze(IF_D(pX,pY,:)),'Linewidth',2);
            xlabel('$z$ in $\lambda$','FontSize',20,'Interpreter','latex');
            ylabel('$I(x=0,y=0,z)$','FontSize',20,'Interpreter','latex')
            set(hs(:),'FontSize',20,'LineWidth',2)
        end
    end

    % HOLOGRAM creation from amplitude and phase
    Amp1 = normalize2D( abs( EDx ));
    Ph1  = angle(EDx);
    Amp2 = normalize2D( abs( EDy ));
    Ph2  = angle(EDy);

    % hologram generator following the Arizon procedure and 
    % following the exp. SLM response (it is the same function as in the lab.)
    [H1,H2] = holoGen(Amp1,Amp2,Ph1,Ph2);

else
% just load the hologram 
    H1 = imread('Hologram1.bmp');
    H2 = imread('Hologram2.bmp');
end
        
if verbose > 0
    figure;
    ax=subplot(1,2,1);
    imagesc(H1( SLM_ROI , SLM_ROI ),[0 255]); cmap('Gray');
    title 'Hologram 1' ; axis square ;
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
            
    ax=subplot(1,2,2);
    imagesc(H2( SLM_ROI , SLM_ROI ),[0 255]); cmap('Gray');
    title 'Hologram 2' ; axis square ;
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
end

%% Arizon procedure Simulation (with or without noise)

if NOISE_flag == 0 % unNoisy case  (ideal SLM calibration)
    wNoise_Amp = [0 0];
    wNoise_Ph  = [0 0];
end

% Getting the complex transmittance of the SLMs:
%  - applying white noise to the response and
%  - filtering in the Fourier Domain [Arizon Procedure Simulation].
[T_SLM1,T_SLM2] = holo_simulator(H1,H2,wNoise_Amp,wNoise_Ph);

Norm = max( max(abs(T_SLM1(:))), max(abs(T_SLM2(:))) );
T_SLM1 = T_SLM1/Norm;
T_SLM2 = T_SLM2/Norm;

switch NOISE_flag 
  case 0   % unNoisy case  (ideal polarizers and mirrors)
    Arm1 = T_SLM1 ;
    Arm2 = T_SLM2 ;
    
  otherwise % Noisy case   (Perlin noise in both Amplitude and Phase)      
    Arm1 = T_SLM1 .* exp( 1i*PerlinNoise(size(T_SLM1),pNoise_Ph) ) ...
               + PerlinNoise(size(T_SLM1),pNoise_Amp);

    Arm2 = T_SLM2 .* exp( 1i*PerlinNoise(size(T_SLM2),pNoise_Ph) ) ...
               + PerlinNoise(size(T_SLM2),pNoise_Amp);
end

% illuminating the system
E0_x = Arm1 .* g .* mask;
E0_y = Arm2 .* g .* mask;

E0x = E0_x( SLM_ROI , SLM_ROI ) ;
E0y = E0_y( SLM_ROI , SLM_ROI ) ;
[f1,f2]=f1f2(E0x, E0y, NA);

if verbose == 0
    figure; 
    if strcmp(ill_FLAG , 'Scalar')
        subplot(1,2,1);
        imagesc(abs2(E0x),[-1 1]); cmap('RedBlue');
        title('$|Ex_0|^2$','Interpreter','Latex') ;
        axis square ; axis off;
        subplot(1,2,2);
        imagesc(abs2(E0y),[-1 1]); cmap('RedBlue');
        title('$|Ey_0|^2$','Interpreter','Latex') ;
        axis square ; axis off;
    else
        subplot(1,2,1);
        imagesc(abs2(f1),[-1 1]); cmap('RedBlue');
        title('$|E_0 \cdot e_1|^2$','Interpreter','Latex') ;
        axis square
        subplot(1,2,2);
        imagesc(abs2(f2),[-1 1]); cmap('RedBlue');
        title('$|E_0 \cdot e_2|^2$','Interpreter','Latex');
        axis square
    end
else
    FontSize = 18;
    figure;
    ax=subplot(2,2,1);
    imagesc(abs2(E0x),[-1 1]); cmap('RedBlue');
    title('$|E_0 \cdot e_x|^2$','Interpreter','Latex') ; 
    axis square ;
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
    
    ax=subplot(2,2,2);
    imagesc(abs2(E0y),[-1 1]); cmap('RedBlue');
    title('$|E_0 \cdot e_y|^2$','Interpreter','Latex')  ; 
    axis square ;
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
    
    ax=subplot(2,2,3);
    imagesc(abs2(f1),[-1 1]); cmap('RedBlue');
    title('$|E_0 \cdot e_1|^2$','Interpreter','Latex')  ; 
    axis square
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
    
    ax=subplot(2,2,4);
    imagesc(abs2(f2),[-1 1]); cmap('RedBlue');
    title('$|E_0 \cdot e_2|^2$','Interpreter','Latex')  ;
    axis square
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
    
    figure; 
    ax=subplot(3,2,1); 
    imagesc(angle(E0x),[-pi pi]);
    title('$\phi_x$','Interpreter','Latex') 
    ;axis square
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
    
    ax=subplot(3,2,2); 
    imagesc(angle(E0y),[-pi pi]);
    title('$\phi_y$','Interpreter','Latex') ;
    axis square
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
    
    ax=subplot(3,2,[3 4 5 6]); 
    imagesc(mod(angle(E0y)-angle(E0x),2*pi),[-pi pi]);
    title('$\phi_y - \phi_x$','Interpreter','Latex') ;
    axis square; cmap('phase');
    ax.FontSize = FontSize;
    ax.LineWidth = 2;
    ax.XTick = 0:60:180;
    ax.YTick = 0:60:180;
    axis([0 180 0 180])
end



%% We can applay an offset for the Optical Axis

if any(disAlign)
    disAlign = abs(disAlign);
    out = zeros(size(E0_x));
    
    out(1:end-disAlign(1),1:end-disAlign(2)) = E0_x(disAlign(1)+1:end,disAlign(2)+1:end);
    out(end-disAlign(1):end,end-disAlign(2):end) = 0;
    E0_x = out;
    
    out(1:end-disAlign(1),1:end-disAlign(2)) = E0_y(disAlign(1)+1:end,disAlign(2)+1:end);
    out(end-disAlign(1)+1:end,end-disAlign(2)+1:end) = 0;
    E0_y = out;
end


%% Richards-Wolf integration for an aplanatic focusing system

[EFx,EFy,EFz] = RWinFFT(E0_x,E0_y,NA,N,L,Zs,'Cartesian');

IF  = squeeze( abs2(EFx(pm:pM,pm:pM,:)) + ...
               abs2(EFy(pm:pM,pm:pM,:)) + ...
               abs2(EFz(pm:pM,pm:pM,:)) );
M_F = max(sqrt(IF(:)));
IF  = IF/M_F/M_F;
EFx_S = EFx(pm:pM,pm:pM,:) / M_F;
EFy_S = EFy(pm:pM,pm:pM,:) / M_F;
EFz_S = EFz(pm:pM,pm:pM,:) / M_F;
IT    = squeeze( abs2(EFx(pm:pM,pm:pM,:)) + abs2(EFy(pm:pM,pm:pM,:)) ) /M_F/M_F;



%% Ploting all together

if ~doSphAbCorr
    Zdefocus = As;
else
    Zdefocus = 0;
end

% cuts at certain Z-plane acording to the Figures on the paper
[~,p_zM] = min(abs(Zs+14-Zdefocus));  % z=-23
[~,p_zF] = min(abs(Zs-0-Zdefocus));  % z=14
[~,p_zm] = min(abs(Zs-14-Zdefocus));  % z=23 

FontSize = 18;
L0 = 4; % final ROI on the XY-cuts 
figure;
ax=subplot(3,3,1);
hp=surf(XY,XY,squeeze(IT(:,:,p_zM)));axis square
set(hp,'edgecolor','none');shading interp;view(2);
xlabel('$x$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
ylabel('$y$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
title(['$|E_t(r,z_1=' num2str(Zs(p_zM),'%4.0f') '\lambda)|^2$'],'Interpreter','Latex','FontSize',FontSize)
axis([-L0 L0 -L0 L0]) 
ax.FontSize = FontSize;

ax=subplot(3,3,2);
hp=surf(XY,XY,squeeze(IT(:,:,p_zF)));axis square
set(hp,'edgecolor','none');shading interp;view(2);
xlabel('$x$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
ylabel('$y$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
title(['$|E_t(r,z_2=' num2str(Zs(p_zF),'%4.0f') '\lambda)|^2$'] , 'Interpreter','Latex','FontSize',FontSize)
axis([-L0 L0 -L0 L0]) 
ax.FontSize = FontSize;

ax=subplot(3,3,3);
hp=surf(XY,XY,squeeze(IT(:,:,p_zm)));axis square
set(hp,'edgecolor','none');shading interp;view(2);
xlabel('$x$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
ylabel('$y$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
title(['$|E_t(r,z_3=' num2str(Zs(p_zm),'%4.0f') '\lambda)|^2$'] , 'Interpreter','Latex','FontSize',FontSize)
axis([-L0 L0 -L0 L0]) 
ax.FontSize = FontSize;

ax=subplot(3,3,[4 5 6]);
hp=surf(Zs,XY,(squeeze(IT(:,pY,:))).^0.8);axis equal
set(hp,'edgecolor','none');shading interp;view(2);
xlabel('$z$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize);
ylabel('$x$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize);
ax.FontSize = FontSize;

ax=subplot(3,3,[7 8 9]);
% hp=surf(Zs,XY,(squeeze(IF(:,pY,:))).^0.8);axis equal
% set(hp,'edgecolor','none');shading interp;view(2);
% xlabel('$z$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize);
% ylabel('$x$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize);
% ax.FontSize = FontSize;

% ax=subplot(3,3,[10 11 12]);
plot( Zs , squeeze(IF(pX,pY,:)) ,'Color',[0.5 0.2 0.2], 'LineWidth',2 )
hold on;
plot(Zs(p_zM),IF(pX,pY,p_zM),'xb' , ...
     Zs(p_zF),IF(pX,pY,p_zF),'xb' , ...
     Zs(p_zm),IF(pX,pY,p_zm),'xb' );
ax.FontSize = FontSize;
ax.LineWidth = 2;
axis([Dz(1) Dz(2) 0 1.05])
xlabel('$z$ (in $\lambda$)','Interpreter','Latex','FontSize',FontSize)
ylabel('$|E(0,z)|^2$','Interpreter','Latex','FontSize',FontSize)
hl=legend('Axis profile','$z_1$','$z_2$','$z_3$');
hl.Interpreter='Latex';
cmap('hot');


%% Finish
toc



