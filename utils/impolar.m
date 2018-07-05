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


% Draw a polar plot with a radius = rho_max. 
% data is a matrix where the 1st dimension is for the radial coordiante and
% the 2nd is for the angular coordinate

function[h]=impolar(data,rho_max,XlabelSTR,YlabelSTR)


phi = linspace(0,  2*pi  ,size(data,2)+1 );
rho = linspace(0,rho_max , size(data,1) );

data(:,end+1)=data(:,1);

[PHI,RHO] = meshgrid(phi,rho);
[ X , Y ] = pol2cart(PHI,RHO);

h=figure('Color',[0.8 0.8 0.8]);

surf(X,Y,data,'edgecolor','none'),
view(0,90) % if data is NOT smooth view(-0.01,90.05)
axis tight
axis equal
shading interp
grid off
xlabel(XlabelSTR)
ylabel(YlabelSTR)
colorbar

