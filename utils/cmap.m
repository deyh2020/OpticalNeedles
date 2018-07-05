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



function[map]=cmap(label)


switch label
    
    case 'RedBlue'
        
        map = [  [linspace(0,1,32)    ones(1,32)   ] ;   ...
                 [linspace(0,1,32) linspace(1,0,32)] ;   ...
                 [   ones(1,32)    linspace(1,0,32)]   ]';
    
    case 'greenW'
        
        map = [  [linspace(1,0,44)   zeros(1,20)] ;   ...
                 [ones(1,44)  linspace(1,0.5,20)] ;   ...
                 [linspace(1,0,44)   zeros(1,20)]   ]';
    
    
    case 'red'  %'red' with blak background is 'hot'
        
        map = [    ones(1,64)                          ;   ...
                 [linspace(1,0,32) linspace(0,1,32)]   ;   ...
                 [linspace(1,0,32)  zeros(1,32)        ]   ]';
             
    case 'green'
        
        map = [  [linspace(1,0,32) linspace(0,1,32)]   ;   ...
                   ones(1,64)                          ;   ...
                 [linspace(1,0,32)  zeros(1,32)        ]   ]';
    

    case 'blueB'
        
        map = [  [  zeros(1,43)  linspace(0,1,21)]   ;   ...
                 [  zeros(1,21)  linspace(0,1,43)]   ;   ...
                 [linspace(0,1,32)   ones(1,32)  ]   ]';         
             
             
             
             
             
    case 'greenB'
        
        map = [  [  zeros(1,21)  linspace(0,1,43)]   ;   ...
                 [linspace(0,1,32)      ones(1,32)     ]   ;   ...
                 [   zeros(1,43)    linspace(0,1,21)   ]   ]';         

             
             
    case 'greenBY'
        
        map = [  [  zeros(1,21)  linspace(0,1,27)  ones(1,16)]   ;   ...
                 [linspace(0,1,21)      ones(1,43)     ]   ;   ...
                 [   zeros(1,48)    linspace(0,1,16)   ]   ]';         

    case 'phase'
        
        map = [  [   linspace(1,0,32)      linspace(0,1,32)      ]   ;   ...
                 [ zeros(1,32) linspace(0,1,16) linspace(1,0,16) ]   ;   ...        
                 [   linspace(0,1,32)      linspace(1,0,32)      ]   ;   ]';
                   
     
    otherwise
        
        map = colormap(label);
        
        
             
             
end
             
             
% axis equal           
             
colormap(map);

% colorbar;
