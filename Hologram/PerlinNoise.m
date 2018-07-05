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


% PerlinNoiseImage generator of size Size1*Size2
% min(PerlinNoiseImage) = delta(1)
% max(PerlinNoiseImage) = delta(2)
%
%% For testing
%   Size = [768,1024]; 
%   intValue = [-.2 .3];
%
function [PerlinNoiseImage] = PerlinNoise(Size,delta)

PerlinNoiseImage = zeros(Size); % output image

w = max(Size);   % width of each iteration
i = 0;           % iteration counter

while w > 3
    i = i + 1;
    d = interp2(randn(w), i-1, 'splines');
    PerlinNoiseImage = PerlinNoiseImage + i*d(1:Size(1),1:Size(2));
    w = w - ceil(w/2-1);
end
  
PerlinNoiseImage = normalize2D(PerlinNoiseImage);
PerlinNoiseImage = delta(1) + PerlinNoiseImage*(delta(2)-delta(1));

%imagesc(PerlinNoiseImage);colorbar