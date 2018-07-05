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



% return the size (height and width) of the CCD in equivalent microns on
% the focal plane for a convencional imaging system, where IM_scale is an
% image which contains the fringes (horizontal and vertical) coorresponding
% to a certaing element-group of a 1951 USAF test imaged with the system to
% be scaled.
%
% This function returns the CCD's full-size instead the equivalent pixel
% size because the capturing software often change the size of the captured
% image, whereas the CCD-size never changes.

function [CCDh_um,CCDw_um] = XY_scale(IM_scale)
% IM_scale = 'USAF_71.png';

USAF_F = 128; % G7-E1: 128 mm-1 (frequency)
USAF_P = 1/USAF_F*1000; % Periodicity in um

IM_scale = im2double(imread(IM_scale));

% size of the image in pixels
N = size(IM_scale);

% ROIs where there are the frinfes
% Horizontal fringes to height scaling
Xi_H = 475; 
Xf_H = 695;
Yi_H = 400;
Yf_H = N(1);
% Vertical fringes to width scaling
Xi_V = 1;
Xf_V = 410;
Yi_V = 480;
Yf_V = 685;

% setting up the figure to show info
figure;
set(gcf, 'Name','Scaling the focal zone transversal units:')
set(gcf, 'Position',[10 150 450 750])

% showing the image
subplot(2,1,1)
imagesc(IM_scale); hold on
plot( [Xi_V Xf_V],[Yi_V Yi_V],'-w' , [Xi_V Xf_V],[Yf_V Yf_V],'-w')
plot( [Xi_H Xi_H],[Yi_H Yf_H],'-w' , [Xf_H Xf_H],[Yi_H Yf_H],'-w')
axis equal;
axis([1 N(2) 1 N(1)])
title('Test 1951 USAF from R3L1S4P-ThorLabs [G7:e1]' , ...
      'FontSize',18 , 'Interpreter','Latex' );
  
annotation( 'textbox',[0.50 0.57 0.4 0.2] , 'String','Height:' , ...
            'LineStyle','none' , 'Interpreter','Latex' , ...
            'FontSize',18 , 'Color',[1 1 1] )       

annotation( 'textbox',[0.25 0.57 0.4 0.2] , 'String','Width:' , ...
            'LineStyle','none' , 'Interpreter','Latex' , ...
            'FontSize',18,'Color',[1 1 1])       

       
% averaging the profiles inside the ROI
profileH = mean( IM_scale( Yi_H:Yf_H , Xi_H:Xf_H )   ,2);
sigmaH   = std(  IM_scale( Yi_H:Yf_H , Xi_H:Xf_H ) ,0,2);
profileW = mean( IM_scale( Yi_V:Yf_V , Xi_V:Xf_V )   ,1);
sigmaW   = std(  IM_scale( Yi_V:Yf_V , Xi_V:Xf_V ) ,0,1);

profileH = profileH/max(profileH);
profileW = profileW/max(profileW);


%% Horizontal Zones [Width]

% find the half maximum of every edge
zones = [1 130 180 230 285 335 410]; % roughly determined zones
for i=1:length(zones)-1

    [~,aux]  = min( abs( profileW(zones(i):zones(i+1)) - 0.5 ) );
    pFWHM_W(i) = aux + zones(i) - 1;
    
end

% averaged fringes period in pixels
PW     = mean(pFWHM_W(3:end)-pFWHM_W(1:end-2));
scaleW = USAF_P/PW; % um/px

% corresponding CCD's width in focal zone (in um)
CCDw_um = scaleW*max(size(IM_scale));



%% Vertical Zones [Height]

% find the half maximum of every edge
zones = [1 70 122 175 221 279 340]; % roughly determined zones
for i=1:length(zones)-1

    [~,aux]  = min( abs( profileH(zones(i):zones(i+1)) - 0.5 ) );
    pFWHM_H(i) = aux + zones(i) - 1;
    
end

% averaged fringes period in pixels
PH     = mean(pFWHM_H(3:end)-pFWHM_H(1:end-2));
scaleH = USAF_P/PH; % um/px

% corresponding CCD's height in focal zone (in um)
CCDh_um = scaleH*min(size(IM_scale)); 



%% ploting profiles and information

subplot(2,1,2);

% --- HORIZONTAL fringes correspondig to the WIDTH of the CCD ---

% ploting the averaged profile +/- standard deviation
[h1,h2] = plotPMsigma( (1:410)' , (profileW)' , (sigmaW)' ); hold on
h2.LineStyle = '-'; 

% ploting the half maximums employed (from up to down)
h3 = plot( pFWHM_W(1:2:end) , profileW(pFWHM_W(1:2:end)) , ...
           'bx' , 'LineWidth',3 , 'MarkerSize',10 );
       
% ploting the half maximums employed (from down to up)
h4 = plot( pFWHM_W(2:2:end) , profileW(pFWHM_W(2:2:end)) , ...
           'rx' , 'LineWidth',3 , 'MarkerSize',10 );


% --- VERTICAL fringes correspondig to the HEIGHT of the CCD ---

% ploting the averaged profile +/- standard deviation
[h5,h6] = plotPMsigma( (420:420+340)' , profileH , sigmaH ); hold on
h6.LineStyle = '-'; 

% ploting the half maximums employed (from up to down)
h7 = plot( pFWHM_H(1:2:end)+420 , profileH(pFWHM_H(1:2:end)) , ...
           'bx' , 'LineWidth',3 , 'MarkerSize',10 );

% ploting the half maximums employed (from down to up)
h8 = plot( pFWHM_H(2:2:end)+420 , profileH(pFWHM_H(2:2:end)) , ...
           'rx' , 'LineWidth',3 , 'MarkerSize',10 );

title 'The 1st element of 7th Group of 1951 USAF test has 128 cycles/mm:'
axis([1 N(1) 0 1.75]);

% annoting some information related with HORIZONTAL fringes [CCD WIDTH]
STR = { ' \quad Width:' , ...
        ['  $\hat{P}_W=' num2str(PW,'%1.2f') 'px$']     , ...
        ['  $\Delta_W=' num2str(scaleW*1E3,2) 'nm/px$'] , ...
        ['  $W_{CCD}=' num2str(CCDw_um,2) '\mu m$' ]    } ;
    
annotation( 'textbox',[0.20 0.24 0.4 0.2] , 'String',STR , ...
            'LineStyle','none', 'Interpreter','Latex' , 'FontSize',18)
 
% annoting some information related with VERTICAL fringes [CCD HEIGHT]
STR = { ' \quad Height:' , ...
        ['$\hat{P}_H=' num2str(PH,'%1.2f') 'px$']     , ...
        ['$\Delta_H=' num2str(scaleH*1E3,2) 'nm/px$'] , ...
        ['$H_{CCD}=' num2str(CCDh_um,2) '\mu m$']     } ;
    
annotation( 'textbox',[0.60 0.24 0.4 0.2] , 'String',STR , ...
            'LineStyle','none' , 'Interpreter','Latex' , 'FontSize',18 )






