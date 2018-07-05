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


% Filtrates the frequencies in the Fourier domain
% Type = 1 : Low-pass ; Type = 2 : High pass

function[filtered]=FourierFilter(image,radius,type)

N = size(image);

[x,y]   = meshgrid(linspace(-N(2)/2,N(2)/2,N(2)),linspace(-N(1)/2,N(1)/2,N(1)));
if     type==1
    filter = ( x.^2+y.^2 <= radius.^2 ).*1;
elseif type==2
    filter = ( x.^2+y.^2 >= radius.^2 ).*1;
end
IM  = fftWELL(image);
IMf = IM.*filter;
filtered = abs(ifftWELL(IMf));