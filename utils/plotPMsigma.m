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



% plots the Y values in the X position with a grayed area sourronding the
% curve indicating the S standard deviation. X, Y, S have to be of the same
% size and must accomplishe size(X)=size(Y)=size(S)=[N 1].
%
% return the handle correspondig to the grayed area (hSigma) and to the
% curve (hY).

function [hSigma,hY] = plotPMsigma(X,Y,S)

% often the values are very sharp and do not look fine, so we smooth it
y1 = smooth(Y-S,4);
y2 = smooth(Y+S,4);
% whereas we want to keep the central values as in the input
y1(ceil(end/2)) = Y(ceil(end/2))-S(ceil(end/2));
y2(ceil(end/2)) = Y(ceil(end/2))+S(ceil(end/2));

% creating the sourranding area
x = [X'  , fliplr(X') ]; % create a continuous x value array for plotting
y = [y1' , fliplr(y2')]; % create y values for out and then back

% ploting the sigma area
hSigma = fill( x,y,'b' ,'FaceColor',[0.7 0.7 0.7],'LineStyle','none'); 
hold on

% ploting the curve
hY = plot( X,Y,'k.' , 'LineWidth',2 , 'MarkerSize',8);

axis([min(X) max(X) 0 max(Y+S)*1.05]);

hold off




