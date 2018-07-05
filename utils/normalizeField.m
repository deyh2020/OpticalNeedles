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


function [E1_out,E2_out] = normalizeField(E1_in,E2_in)

E1_out = E1_in/max([max(E1_in(:)) max(E2_in(:))]);
E2_out = E2_in/max([max(E1_in(:)) max(E2_in(:))]);