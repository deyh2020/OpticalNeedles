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


% Gets an image and returns the angular averaged value Vavg that has a
% radial distance Rho from the [centerX centerY] and a standard deviation
% Savg associated to the dispersion of the averaged Navg values.

function [Rho,Vavg,Savg,Navg] = angularAverage (image,centerX,centerY)
    
N = size(image);

valueVSrho = zeros(N(1)*N(2),2);

auxSC = 1E8; % to get decimals enough in float2integer
count = 1;

for i = 1:N(1)
     for j = 1:N(2)
         
         valueVSrho(count,1) = round(  ...               % integer
            sqrt( (j-centerY).^2+(i-centerX).^2 )  ...  % rho
                                             *auxSC );    % with 8 decimals
         valueVSrho(count,2) = image(i,j);        % Value of certain pixel

         count = count+1;
     end
end

valueVSrho(:,:)=sortrows(valueVSrho(:,:)); %Sorting for rhos

data  = nan(N(1)*N(2),4);
count = 1;
rhold = 1542345; % any number to initializate (1 or 0 is not recomended)

for rho=valueVSrho(:,1)' % Scanning all rhos
    if rho~=rhold        % to avoid any rho already computered

        pp    = (valueVSrho(:,1)==rho); % where are same rho's
        first = find(pp,1,'first');
        last  = find(pp,1,'last');
        Vmean = mean(valueVSrho(first:last,2)); % meaning the values
        Vstd  = std(valueVSrho(first:last,2));  % standard deviation

        data(count,1) = rho;
        data(count,2) = Vmean;
        data(count,3) = Vstd;
        data(count,4) = last-first+1;

        count=count+1;
        rhold=rho;
    end
end

Nd = count-1; % 

% from stricly positive Rho to a symetric coord Rho=[-max(rho) max(rho)]
Rho  = [-data(Nd:-1:2,1);data(1:Nd,1)]/auxSC; % integer2float with decimals
Vavg = [ data(Nd:-1:2,2);data(1:Nd,2)];  % Angular averaged values
Savg = [ data(Nd:-1:2,3);data(1:Nd,3)];  % Standard Deviation dispersion
Navg = [ data(Nd:-1:2,4);data(1:Nd,4)];  % Number of values employed


