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


%   ANTIMOD Avoid jumps in multivalued functions
%
% [out,pFirst]=ANTIMOD(in,P) avoids jumps in multivalued functions as 
% angular/phase functions 'out' is a continous array based on 'in' avoiding
% jumps of 'P' points (or multiples of 'P'). 'pFirst' is the index of the 
% first avoided jump and 'N' is the number of jumps corrected.

function[out,pFirst]=antimod(in,P)

N      = 0;
pFirst = 0; 

for i = 2:length(in)
    
    Delta = in(i-1)-in(i);
        
    if abs(Delta) > P*0.5
        in(i) = in(i) + P*sign(Delta)*round(Delta/P);
        if N==0, pFirst = i; end
        N = N+1;
    end
    
end

out = in;

disp(['Antimod has corrected ' num2str(N) ' points.'])