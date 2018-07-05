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



% function[]=delipse(S0,S1,S2,S3)

% Entren les componets complexes (x,y) del camp. Està pensat per 16x16.
% Selecciona els 8x8 de la part central.
% Np és el nombre de punts en que dibuixa cada el·lipse.

function[imOUT,hs,hhh]=delipse(S0,S1,S2,N)

Np  = round(5000/N/N);
Nim = max(size(S0));
s0  = zeros(N,N);
s1  = zeros(N,N);
s2  = zeros(N,N);
%s3  = zeros(N,N);

for i=1:N
    for j=1:N
        int_i   = floor((i-1)*Nim/N+1):floor(i*Nim/N);
        int_j   = floor((j-1)*Nim/N+1):floor(j*Nim/N);
        
        s0(i,j) = mean(mean( S0(int_i,int_j) ));%S0(round((i-0.4)*Nim/N),round((j-0.4)*Nim/N));%
        s1(i,j) = mean(mean( S1(int_i,int_j) ));%0.6*S1(round((i-0.4)*Nim/N),round((j-0.4)*Nim/N));%
        s2(i,j) = mean(mean( S2(int_i,int_j) ));%0.6*S2(round((i-0.4)*Nim/N),round((j-0.4)*Nim/N));%
        %s3(i,j) = mean(mean( S3(int_i,int_j) ));%S3(round((i-0.4)*Nim/N),round((j-0.4)*Nim/N));%
    end
end
     
% L=sqrt(s1.^2+s2.^2);
A0=(s0-s1);%sqrt((s0+L)/2);
B0=(s0+s1);%sqrt((s0-L)/2);
%theta=atan2(s2,s1)/2;
M=max(max(max(A0)),max(max(B0)));
A=A0/max(max(M))/2.25;
B=B0/max(max(M))/2.25;
Cph=-s2./(sqrt(s0.^2-s1.^2-s2.^2));
%Sph=s3./sqrt(s0.^2-s1.^2-s2.^2);
ph=acos(Cph);
%phs=asin(Sph);



hf = figure;
hi = imagesc(linspace(0.5,N+0.5,Nim),linspace(0.5,N+0.5,Nim),S0);
hs = {hi};
hhh = {[0 0]};
cmap('gray');
hold on
x=zeros(1,Np);
y=zeros(1,Np);
for ii=1:N;
    for jj=1:N;
        for t1=0:Np;
            t=t1*2*pi/Np;
%             x=ii+A(ii,jj)*cos(t)*cos(theta(ii,jj))-B(ii,jj)*sin(t)*sin(theta(ii,jj));
%             y=jj+A(ii,jj)*cos(t)*sin(theta(ii,jj))-B(ii,jj)*sin(t)*cos(theta(ii,jj));
            y(t1+1)=real(ii+A(ii,jj).*cos(t));
            x(t1+1)=real(jj+B(ii,jj).*cos(t+ph(ii,jj)));
        end
        hp = plot(x,y,'-g','LineWidth',2);
        hs{end+1} = hp;
        hhh{end+1} = [ii jj];
    end
end
axis([0 N+1 0 N+1]);
axis off;
axis square;
hold off;
F=getframe(hf);
[imOUT,Map]=frame2im(F);

