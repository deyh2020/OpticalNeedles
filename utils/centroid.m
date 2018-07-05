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


% Search the centroid of an image
function[row,column]=centroid(matrix)

MAX = max(matrix(:));
matrix(matrix<0.5*MAX) = 0;
[rc,cc] = meshgrid(1:size(matrix,1),1:size(matrix,2));
Mt      = sum(matrix(:));
column  = round(sum(sum(  matrix.*cc.'  )) / Mt);
row     = round(sum(sum(  matrix.*rc.'  )) / Mt);
