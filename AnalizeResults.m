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




close all;
clear variables
tic

%% PATHs and info LOADING 

addpath(genpath('Holograms'))
addpath(genpath('utils'))

currentPATH = cd;
currPathSep = strsplit(currentPATH,'/');
pRootPATH   = find( strcmp(currPathSep,'BeamsModelationRepo') );
SystemPATH  = '';

for i=1:pRootPATH
    SystemPATH = [SystemPATH '/' currPathSep{i}];
end

SystemPATH  = [SystemPATH '/DigiHolos2LaserBeamModelation'];

fileID   = fopen([SystemPATH '/+scripts/design_kinds.txt']);
Bkinds{1} = '';
Btxt{1}   = '';
i=1;
while ischar(Bkinds{end})
    string = fgets(fileID);
    Bkinds{i} = string(1:end-2);
    Blbl(i) = str2double(string(1:find(string==':')-1));
    Btxt{i} = string(find(string==':')+2:end);
    i=i+1;
end
fclose(fileID);



%% Loading results files

dirBAUL = '200317_51';
fol     = dir([dirBAUL '*']);
Folders = fol([fol.isdir]);

for iFolder = 1 : length(Folders)
Folder = Folders(iFolder);
    
load([Folder.name '/sim_N8_m8_esc.mat']);
Zs_Theo = Zs;
beam_label = str2double(Folder.name(end-1:end));
    beam_str   = Btxt{beam_label};
    pN = 12;
    pm = 17;
    N_beam     = str2double(beam_str(pN+2:pN+3));
    m_beam     = str2double(beam_str(pm+2:pm+3));
    if beam_label > 66, As = str2double(beam_str(end-1:end)); else As=0; end
    disp(['N=' num2str(N_beam) ' ; m=' num2str(m_beam) ' ; As=' num2str(As)])

    

%% Parameters and variables

centroid_FLAG  = 3; % 1:centroid ; 0:peak ; 2:peak&centering.m ; 3:certain value
centroid_CHECK = 0; % draw ALL images to check the center assignation

XYunits   = '\lambda';  % chosse between [ 'nm' , 'um' , '\lambda' ]
Zunits    = '\lambda';  % chosse between [ 'nm' , 'um' , '\lambda' ]
lambda    = 594;        % nm
FILEunits = 1;          % 1 for nm ; 1000 for un (fileName contains z-position)

% equivalent CCD size in focal area [um]
[CCDscaleH,CCDsclaeW] = XY_scale('USAF_71.png'); 

L = 70;  % half side of ZOI in PIXELS

%% Loading files and centering the beam
%  ims(Nx,Ny,n) is the 3D matrix with all the insets images
%  Znames(n) is Z position acording to the name of the file
%  center(2,n) is the maxima position on the images
tic
display(' -- Loading files: -- ')

% Prepare the new file.
vidObj = VideoWriter([Folder.name '.avi']);
vidObj.FrameRate = 5;
open(vidObj);

% Create an animation.
figure1 = figure;
imagesc(zeros(2*L+1,2*L+1));
colormap gray;
axis equal
axis tight
set(gca,'nextplot','replacechildren');

currentdir = cd;
cd(Folder.name);
files = dir('*.png');
Zs    = max(size(files));  % number of files = number of z-cuts
Nim   = size(imread(files(1).name));
CCDh  = min(Nim);
scale = CCDscaleH/CCDh; % scale in um/px

% initializing
Zname   = zeros(Zs,length(files(1).name));
Zcurr   = zeros(1,Zs); % doubles in nm
centers = zeros(2,Zs); % center of each plane
ims     = zeros(L*2+1,L*2+1,Zs);
imsT    = zeros([size(imread(files(1).name)),Zs]);

[x,y]   = meshgrid( linspace(1,Nim(2),Nim(2)) , ...
                    linspace(1,Nim(1),Nim(1))  );

i=1; % initializing the counter
for file=files'

    im         = im2double(imread(file.name));
    imf        = FourierFilter(im,60,1);
    imf        = im/max(im(:));
    %im(:,1:200)= 0;
    Zname(i,:) = file.name(1:end);
    Zcurr(i)   = str2double(file.name(1:8))*FILEunits; % in nm   

    
    if     centroid_FLAG == 1  

        M = im;
        %M(:,1:300)=0;
        [Cx,Cy]=centroid(M); 
        
    elseif centroid_FLAG == 3

        [Cx,Cy]=centroid(im); 

        centering;  % experimental dir should have a centering.m file
        centers(:,i)=[Cx;Cy];
        
    elseif centroid_FLAG ~= 1
        
        [~,centers(1,i)] = max(max(im(321-L:321+L , 512-L:512+L  )));
        [~,centers(2,i)] = max(max(im(321-L:321+L , 512-L:512+L  ),[],2));
        
        if centroid_FLAG == 2
            centering;
        end
        if any([centers(1,i)+L>size(im,2) ...
                centers(2,i)+L>size(im,1) ])
            centers(2,i) = 354;
            centers(1,i) = 475;
            disp(num2str(i))
        end
        Cx = centers(2,i)+321-L;
        Cy = centers(1,i)+512-L;
    end
    
    ims(:,:,i) = im(Cx-L:Cx+L , Cy-L:Cy+L  );
    
    %  Drawing a cercle sourronding the center (just for visualizing)
    imsV = imf/max(imf(:));            
    imsV( Cx-5 , Cy   ) = 1;
    imsV( Cx+5 , Cy   ) = 1;
    imsV( Cx   , Cy+5 ) = 1;
    imsV( Cx   , Cy-5 ) = 1;
    imsV( Cx-3 , Cy-3 ) = 1;
    imsV( Cx+3 , Cy-3 ) = 1;
    imsV( Cx-3 , Cy+3 ) = 1;
    imsV( Cx+3 , Cy+3 ) = 1;
    
    % relative z-pos
    Zplot = (Zcurr(i)-Zcurr(1))/lambda;             
    figure(figure1);
    imagesc(normalize2D(imsV(Cx-L:Cx+L , Cy-L:Cy+L  )));
    colormap gray;axis equal

    if i==1
        % filename annotation
        annot1 = annotation(figure1,'textbox',...
                            [0.3 0.14 0.7 0.08],...
                            'Color',[1 1 1],...
                            'String',['file(' num2str(i) '): ' Zname(i,:)],...
                            'LineStyle','none',...
                            'FontWeight','bold',...
                            'FontSize',24,...
                            'FitBoxToText','off');

        % z-pos annotation
        annot2 = annotation(figure1,'textbox',...
                            [0.4 0.8 0.7 0.08],...
                            'Color',[1 1 1],...
                            'String',['z = ' num2str(Zplot) ' \lambda '],...
                            'LineStyle','none',...
                            'FontWeight','bold',...
                            'FontSize',24,...
                            'FitBoxToText','off');
    else
        % updating the annotations
        annot1.String = ['file(' num2str(i) '): ' Zname(i,:)];
        annot2.String = ['z = ' num2str(Zplot) ' \lambda '];
    end
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
    
    
    if centroid_CHECK == 1 && mod(i-1,25)==0
        % we just check some images to debug
        imsT(:,:,i)= imf;            
        imsT( Cx-5 , Cy   ,i) = 1;
        imsT( Cx+5 , Cy   ,i) = 1;
        imsT( Cx   , Cy+5 ,i) = 1;
        imsT( Cx   , Cy-5 ,i) = 1;
        imsT( Cx-3 , Cy-3 ,i) = 1;
        imsT( Cx+3 , Cy-3 ,i) = 1;
        imsT( Cx-3 , Cy+3 ,i) = 1;
        imsT( Cx+3 , Cy+3 ,i) = 1;
        figure;
        imagesc(imsT(:,:,i));
        title(['file(' num2str(i) '): ' Zname(i,:)]);
    end
    i=i+1;
end
cd(currentdir)

close(vidObj)
close(figure1)

toc;display(' -------------------- ');display(' ')


%% Size and units

N=L*2+1;

switch XYunits
    case 'nm'
        SC=scale*1000;
    case 'um'
        SC=scale;
    case '\lambda'
        SC=scale*1000/lambda;
    otherwise
        error('chosse between [ nm , um , \lambda ] for the XYunit')
end

switch Zunits
    case 'nm'
        Zcurr=Zcurr;
    case 'um'
        Zcurr=Zcurr/1000;
    case '\lambda'
        Zcurr=Zcurr/lambda;
    otherwise
        error('chosse between [ nm , um , \lambda ] for the Zunit')
end

%initializes
Irho    = zeros(N*N,2,Zs); % Irho(:,rho:I,z)
Irho2   = zeros(N*N,2,Zs); % Irho(:,rho:I,z)
if L==100,Nnew=3737;else Nnew=1;end
data    = zeros(Nnew,4,Zs);

%% From images to array ; I vs rho ; doing meanings for all directions
tic
disp(' -- Procesing data: -- ')

for z=1:Zs

    [DATA(:,1,z),DATA(:,2,z),DATA(:,3,z),~]=angularAverage(ims(:,:,z),L+1,L+1);
    
    if mod(z,round(Zs/10) )==0
        disp([num2str(ceil(z/Zs*100)) '% of angular average -> ' ...
              num2str(round(toc)) 's elapsed'])
    end   
end

Nd2 = size(DATA);


toc;display(' --------------------- ');display(' ')


%% Plot the info (disabled)

Zstep=ceil(Zs/10);
for z=1:Zstep:Zs
    
    interval = find(DATA(:,1,z),1,'first') : find(DATA(:,1,z),1,'last');
    Ik  = DATA(interval,2,z)-min(DATA(interval,2,z));
    Ikk = Ik/max(Ik);
    figure
    plot(DATA(interval,1,z)*SC,Ikk)
    title(['Z = ' num2str(Zcurr(z)) ' ' Zunits])
    xlabel(['\rho (' XYunits ')'])
end

    
%% Waist evaluation

thW=1/2;%exp(-2);

W = zeros(1,Zs)*20;
for z=1:Zs
    [~,p] = min(abs( DATA(1:ceil(Nd2(1)/2),2,z)/max(DATA(1:ceil(Nd2(1)/2),2,z))-thW ));
    W(z)  = -DATA(p,1,z)*SC;
end

% centering Z in the minimum waist
[W0,pZ0] = min(W);
Z0       = Zcurr(pZ0);
Z        = Zcurr-Z0;

W_m = mode(W);

figure;hold on
plot(Z,  W, '.', 'MarkerSize', 12);
plot(Z, -W, '.', 'MarkerSize', 12);

xlabel(['$z$ ($'  Zunits '$)'] , 'Interpreter','LaTeX')
ylabel(['$W$ ($' XYunits '$)'] , 'Interpreter','LaTeX')

title(['$FWHM = ' num2str(W0*2,3) XYunits '$ ;  ' ...
       'mode$\{FWHM\} = ' num2str(W_m*2,3) XYunits '$'] , 'Interpreter','LaTeX');
axis([min(Z) max(Z) -max(W)*.25 max(W)*.25])
set(gca, 'FontSize',18 , 'FontWeight','bold' , 'TickLabelInterpreter','LaTeX', ...
         'LineWidth',2);


%% ZR plane formation


Rsize   = round(L*sqrt(2)); %half side (resolution) in pixels for rho coordinate
R       = zeros(1,Rsize*2);
ZRplane = zeros(Rsize*2,Zs);
MAXrho  = max(max(DATA(:,1,:)));
Drho    = max(MAXrho)/Rsize;

Ener = zeros(   1  ,Zs);
I    = zeros(Nd2(1),Zs);
Tin  = zeros(   1  ,Zs);
aux0 = 0 ; aux1 = 0 ; aux2 = 0 ; aux3 = 0 ; aux4 = 0;
for z=1:Zs
    
     Ener(z) = trapz(DATA(:,1,z),DATA(:,2,z));
     I1(:,z) = squeeze(DATA(:,2,z))/Ener(z);
     I2(:,z) = squeeze(DATA(:,2,z))/max(squeeze(DATA(:,2,z)));
%     if     Zcurr(z) < 19556.25
%         I(:,z) = squeeze(DATA(:,2,z));
%         aux0 = aux0+1;
%     elseif Zcurr(z) < 19563.75
%         I(:,z) = squeeze(DATA(:,2,z))*1.1;
%         aux1 = aux1+1;
%     elseif Zcurr(z) < 19568.75
%         I(:,z) = squeeze(DATA(:,2,z))*1.2;
%         aux2 = aux2+1;
%     elseif Zcurr(z) < 19573.75
%         I(:,z) = squeeze(DATA(:,2,z))*1.1;
%         aux3 = aux3+1;
%     else
%         I(:,z) = squeeze(DATA(:,2,z));
%         aux4 = aux4+1;
%     end
end
save([Folder.name '/ZRplane.mat'],'I1','beam_label')

I1 = I1/max(I1(:));
I3 = squeeze(DATA(:,2,:));
I3(:,212:241) = I3(:,212:241)*2.5;
I3 = I3/max(I3(:));


%I = imgaussfilt(I,1.4);
I=I1;
Ener2 = zeros(1,Zs);
for z=1:Zs
    Ener2(z) = trapz(DATA(:,1,z),I(:,z));
end
figure;
plot(Ener2)

% centering Z in the maximum in-axis
[~,pZ0]  = max(I(ceil(Nd2(1)/2),:));
Z0       = Zcurr(pZ0);
Z        = Zcurr-Z0;


for z=1:Zs
    
    ZRplane(Rsize,z) = I(ceil(Nd2(1)/2),z); % optical axis
    R(Rsize)         = 0;
    
    for r=1:Rsize-1                % off-axis
        pp    = (DATA(:,1,z)>=r*Drho & DATA(:,1,z)<r*Drho +Drho);
        first = find(pp,1,'first');
        last  = find(pp,1,'last');
        N(r)  = last-first; 
        m     = mean(I(first:last,z));
        
        ZRplane(Rsize+r,z) = m;
        ZRplane(Rsize-r,z) = m;
        R(Rsize-r)         = - (r*Drho+Drho/2)*SC;
        R(Rsize+r)         =   (r*Drho+Drho/2)*SC;
    end
    
end
R(end)=R(end-1)+Drho/2*SC;


%ZRplane = imgaussfilt(ZRplane,1.4);

figure;
%contour(Z,R,normalize2D(ZRplane),[0.5 0.5])%,'LineStyle','none');shading interp
imagesc(Z,R,ZRplane)
xlabel(['$z$ ($' Zunits '$)'] , 'Interpreter','LaTeX')
ylabel(['$r$ ($' XYunits '$)'] , 'Interpreter','LaTeX')
title( Bkinds{beam_label} , ...
           'Interpreter','LaTeX');
set(gca, 'FontSize',30 , 'FontWeight','bold' , 'TickLabelInterpreter','LaTeX');
cmap('hot');
axis equal
axis([min(Z) max(Z) min(R) max(R) 0 1]);view(2)
end


%% Show certain XYplane z=-23,14,23 lambdas

close all

[~,p1] = min(abs( -23 - Z) );
[~,p2] = min(abs( 14 - Z) );
[~,p3] = min(abs( 23 - Z) );

SCbar = floor(1/SC/lambda*1E3); %pix/XYunits

im1 = normalize2D( ims(25:end-25,25:end-25,p1) );
im2 = normalize2D( ims(25:end-25,25:end-25,p2) );
im3 = normalize2D( ims(25:end-25,25:end-25,p3) );
for iZ = 1:3
    pS = [p1 p2 p3];
    [~,pZ] = min( abs( Zs_Theo - Z(pS(iZ)) ));
    [Ravg(:,iZ), Iavg(:,iZ), Savg(:,iZ), ~] = angularAverage( ...
                                   normalize2D(ims(25:end-25,25:end-25,pS(iZ))), ...
                                   size(im1,1)/2+1,size(im1,2)/2+1);
end

Ravg = Ravg*SC;

im1(end-9:end-7,8:8+SCbar-1) = 1;
im2(end-9:end-7,8:8+SCbar-1) = 1;
im3(end-9:end-7,8:8+SCbar-1) = 1;

% figure;
% imagesc(im1);axis square;cmap('hot');axis off
% figure;
% imagesc(im2);axis square;cmap('hot');axis off
% figure;
% imagesc(im3);axis square;cmap('hot');axis off

imSSS(:,:,1) = im1;
imSSS(:,:,2) = im2;
imSSS(:,:,3) = im3;


plotheight=20;
plotwidth=16;
subplotsx=3;
subplotsy=3;   
leftedge=1.5;
rightedge=1.5;   
topedge=2;
bottomedge=1;
spacex=0.5;
spacey=0.6;
fontsize=15;    
 
sub_pos = subplot_pos(plotwidth, plotheight, ...
                      leftedge, rightedge, bottomedge, topedge, ...
                      subplotsx, subplotsy, spacex, spacey);
 
sub_pos{1,3}(2) = sub_pos{1,3}(2)+0.06;
sub_pos{2,3}(2) = sub_pos{2,3}(2)+0.06;
sub_pos{3,3}(2) = sub_pos{3,3}(2)+0.06;

%sub_pos{1,1}(3) = sub_pos{1,1}(3)*3.3;
%sub_pos{1,1}(2) = sub_pos{1,1}(2)-0.05;
sub_pos{1,1}(3) = sub_pos{1,1}(3)*3.3;
sub_pos{1,1}(4) = sub_pos{1,1}(4)*1.1;

%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf, 'Name','Figure4')
set(gcf, 'Color',[1 1 1]);

STR_Zs = { ['$z=' num2str(round(Z(p1))) '\lambda$'] , ...
           ['$z=' num2str(round(Z(p2))) '\lambda$'] , ...
           ['$z=' num2str(round(Z(p3))) '\lambda$'] }  ;
STR_pos={'left','right','right'};

%loop to create axes
count = 1;
for i=1:subplotsx;
    for ii=1:subplotsy;
     
        ax = axes('position', sub_pos{i,ii}, 'XGrid','off' , 'XMinorGrid','off' , ...
                  'FontSize',fontsize , 'Box','on' , 'Layer','top');

        if ii == 2 
            imagesc(imSSS(:,:,i));
            axis tight
            axis off
            axis square
            annotation(f, 'textbox',sub_pos{i,ii} , 'String',[STR_Zs{i} ' '], ...
                          'Interpreter','Latex' , 'LineStyle','none' , ...
                          'FontSize',24 , 'Color',[1 1 1] , 'FontWeight','bold' , ...
                          'HorizontalAlignment','right' , 'VerticalAlignment','bottom')
        end  % ii=2

        if ii == 3
            [~,pZ] = min( abs( Zs_Theo - Z(pS(i)) ));

            [h1,h3,hl] = plotPMsigma(Ravg(:,i), ...
                                Iavg(:,i)-0.1)/max(Iavg(:,i)-0.1),Savg(:,i)/max(Iavg(:,i));
            hold on
            h4 = plot(XY,squeeze(IF_Dt(pX,:,pZ))/max(squeeze(IF_Dt(pX,:,pZ))), ...
                                 'b-','Linewidth',3);hold on
            h5 = plot(XY,squeeze( IF_D(pX,:,pZ))/max(squeeze(IF_Dt(pX,:,pZ))), ...
                                 'r--','Linewidth',2);

            axis([-2.5 2.5 -0.1 max(squeeze( IF_D(pX,:,pZ))/max(squeeze(IF_Dt(pX,:,pZ))))*1.05]);
            hold off

            annotation(f, 'textbox',sub_pos{i,ii} , 'String',STR_Zs{i} , ...
                          'Interpreter','Latex' , 'LineStyle','none' , ...
                          'FontSize',20 , 'Color',[0 0 0] , 'FontWeight','bold' , ...
                          'HorizontalAlignment',STR_pos{i} , 'VerticalAlignment','top')
            
            ax.FontSize=20;
            ax.YTick=([0 0.5 1]);
            ax.YTickLabels=({'0','0.5','1'});
            ax.XTick=(-2:2);
            ax.XTickLabels=({'2','1','0','1','2'});
            ax.XLabel.String = '$r$ in $\lambda$';
            ax.XLabel.Interpreter = 'LaTeX';
            ax.YLabel.String = 'Irradiance';
            ax.YLabel.Interpreter = 'LaTeX';
            ax.LineWidth=2;
            if i==1
                legend off
            end
            if i==2
                ax.YTickLabels=({' ',' ',' '});
                ax.YLabel.String = '';
                hl=legend([h3 h1 h4 h5],...
                    'Experimental','$\pm\sigma$','Theo. Trans.','Theo. Total');
                hl.Interpreter='LaTeX';
                hl.Position=[0.25 0.82 0.2 0.2];
                legend off
            end
            if i==3
                ax.YAxisLocation='right';
                legend off
            end
                

            hold off
        end  % ii=3

        if ii == 1
            if i==1
                surf(Z,R,ZRplane*10);cmap('red');%view(2);
                shading interp;axis equal;view(2);
                xlabel('$z$ in $\lambda$','FontSize',24,'Interpreter','latex');
                ylabel('$r$ in $\lambda$','FontSize',24,'Interpreter','latex'); 
                set(gca, 'FontSize',24 , 'FontWeight','bold' , 'TickLabelInterpreter','LaTeX');
                ax.FontSize=24;
                ax.YTick=([-10 -5 0 5 10]);
                ax.YTickLabels=({'10' '5' '0','5','10'});
                ax.XTick=(-40:20:40);
                %ax.XTickLabels=({'2','1','0','1','2'});
                ax.YLabel.Interpreter = 'LaTeX';
                ax.XLabel.Interpreter = 'LaTeX';
            else
                axis off
            end
        end  % ii=1

    end  % loop for ii
end  % loop for i
cmap('hot');



%% 

Rs = DATA(:,1,z)*SC; 
I_filtered = imgaussfilt(squeeze(I1),4);
figure;
surf(squeeze(Z),squeeze(Rs),normalize2D(squeeze(I1))*10,'LineStyle','none');
view(2);shading interp;cmap('hot');%
axis equal;
xlabel(['$z$ ($' Zunits '$)'] , 'Interpreter','LaTeX')
ylabel(['$\rho$ ($' XYunits '$)'] , 'Interpreter','LaTeX')
%title(['W_0 = ' num2str(W0) ' ' XYunits '  ;  NA = ' num2str(NA) ...
%    ' ; file(z_0) = ' Zname(pZ0,:)]);
%axis([-30 30 -18 18])
title(['$N = ' num2str(N_beam) '$ ; $m = ' num2str(m_beam) '$'] , 'Interpreter','LaTeX');
set(gca, 'FontSize',18 , 'FontWeight','bold' , 'TickLabelInterpreter','LaTeX');


% % figure;plot(Z,I(end/2,:))
% % xlabel(['z (' Zunits ')'])
% % ylabel 'I(\rho=0,z)'




% 
% Z_ws(:,BeamKind)  = Z;
% Rs_ws(:,BeamKind) = Rs;
% I_ws(:,:,BeamKind)= I;
% W_ws(:,BeamKind)  = W;

%%

figure
plot( Z,I(ceil(Nd2(1)/2),:) , 'LineWidth',3)
xlabel(['$z$ ($' Zunits '$)'] , 'Interpreter','LaTeX')
ylabel(['$\rho$ ($' XYunits '$)'] , 'Interpreter','LaTeX')
set(gca, 'FontSize',18 , 'FontWeight','bold' , 'TickLabelInterpreter','LaTeX','LineWidth',2);
axis([min(Z) max(Z) 0 1.05*max(I(ceil(Nd2(1)/2),:))]);


%% image certain file

for i=229:247
figure;
imagesc(imread([dirBAUL '/' files(i).name]))
title(num2str(i))
end
