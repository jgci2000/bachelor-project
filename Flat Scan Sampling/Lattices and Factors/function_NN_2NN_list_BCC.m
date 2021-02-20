function [NN1_zpos1,NN1_zpos2,NN1_zpos3,NN1_zpos4,NN1_zneg1,NN1_zneg2,NN1_zneg3,NN1_zneg4, NN2_xpos,NN2_xneg,NN2_ypos,NN2_yneg,NN2_zpos,NN2_zneg] = function_NN_2NN_list_BCC(L,i,j,k)
%
% LATTICE CONVENTION:
% i = rows , j = columns , k = layers
% top row  =  lowest i value
% left column = lowest j value
% bottom layer = lowest k value
%
% 1NN (BCC)
%

NN1_zpos1 = 666;
NN1_zpos2 = 666;
NN1_zpos3 = 666;
NN1_zpos4 = 666;
NN1_zneg1= 666;
NN1_zneg2 = 666;
NN1_zneg3 = 666;
NN1_zneg4 = 666;

%
% CORE
%
if i~=1 && i~=L && j~=1 && j~=L && k~=1 && k~=2*L %%% WORKS ONLY FOR ODD k?
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end
%

%
% #20, #52, #84
%
if i==1 && j==L && k~=1 && k~=2*L && mod(k,2)==0 
%
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% #29, #61, #93
%
if i==L && j==1 && k~=1 && k~=2*L && mod(k,2)==0 
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% #36, #68, #100
%
if i==1 && j==L && k~=1 && k~=2*L && mod(k,2)==1  
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 ;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L -1; 
    %
end

%
% #45, #77, #109
%
if i==L && j==1 && k~=1 && k~=2*L && mod(k,2)~=0
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2 -1; 
    %
end

%
% Y-BOTTOM LAYER - PLANE CORE, EVEN k #62
%
if i==L && j~=1 && j~=L && k~=1 && k~=2*L && mod(k,2)==0 
 %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% Y-BOTTOM LAYER - PLANE CORE, ODD k #46 (CORE-LIKE)
%
if i==L && j~=1 && j~=L && k~=1 && k~=2*L && mod(k,2)~=0 
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% Y-TOP LAYER - PLANE CORE #34, #50
%
if i==1 && j~=1 && j~=L && j~=1 && k~=2*L
 %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L -1; 
    %
end

%
% X-LEFT LAYER - PLANE CORE #37, #53
%
if i~=1 && i~=L && j==1 && k~=1 && k~=2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2 -1; 
    %
end
%
% X-RIGHT LAYER - PLANE CORE #56 (EVEN), #72 (ODD)
%
if i~=1 && i~=L && j==L && k~=1 && k~=2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% Z-TOP LAYER - PLANE CORE #118
%
if i~=1 && i~=L && j~=1 && j~=L && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2;
    %
end
%

%
% Z-BOTTOM LAYER - PLANE CORE
%
if i~=1 && i~=L && j~=1 && j~=L && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -L -1;
    %
    NN1_zneg1 = 2*L^3 + j+(i-1)*L+(k-1)*L*L -L^2 -L;
    NN1_zneg2 = 2*L^3 + j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = 2*L^3 + j+(i-1)*L+(k-1)*L*L -L^2 -1;
    NN1_zneg4 = 2*L^3 + j+(i-1)*L+(k-1)*L*L -L^2 -L -1;
    %
end



%
% LEFT CORE, TOP LAYER #117, #124
%
if i~=1 && i~=L && j==1 && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% BOTTOM CORE, TOP LAYER #126, #127
%
if i==L && j~=1 && j~=L && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% RIGHT CORE, TOP LAYER #120, #124
%
if i~=1 && i~=L && j==L && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 -L +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2  +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end
%

%
% TOP CORE, TOP  LAYER #114, #115
%
if i==1 && j~=1 && j~=L && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% LEFT CORE, BOTTOM LAYER #5, #9
%
if i~=1 && i~=L && j==1 && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 ;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 +L -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -1; 
    %
end


%
% BOTTOM CORE, BOTTOM LAYER #14, #15 <<<<<<<<<<<<<<<<<<<<<<
%
if i==L && j~=1 && j~=L && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -L -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L -1; 
    %
end

%
% RIGHT CORE, BOTTOM LAYER #8, #12
%
if i~=1 && i~=L && j==L && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -L -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 ;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L -1; 
    %
end

%
% TOP CORE, BOTTOM LAYER #2, #3
%
if i==1 && j~=1 && j~=L && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 ;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L -1; 
    %
end


%
% TOP LEFT, TOP LAYER #113
%
if i==1 && j==1 && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% BOTTOM LEFT, TOP LAYER #125
%
if i==L && j==1 && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% BOTTOM RIGHT, TOP LAYER #128
%
if i==L && j==L && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 -L +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -2*L^2 +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2; 
    %
end

%
% TOP RIGHT, TOP LAYER #116
%
if i==1 && j==L && k==2*L
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 -L +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L -2*L^3 +L^2; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2 ; 
    %
end

%
% BOTTOM RIGHT, BOTTOM LAYER #16
%
if i==L && j==L && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -L -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L -1; 
    %
end


%
% BOTTOM LEFT, BOTTOM LAYER #13
%
if i==L && j==1 && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -1; 
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 ;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 +L -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L^2 -1; 
    %
end

%
% TOP RIGHT, BOTTOM LAYER #4
%
if i==1 && j==L && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L -1;
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +2*L^3;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +2*L^3 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L -1;
    %
end

%
% TOP LEFT, BOTTOM LAYER #1
%
if i==1 && j==1 && k==1
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +2*L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +2*L^2 -1;
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L +2*L^3 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +2*L^3;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +2*L^3 +L -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L +2*L^3 -1;
    %
end

%
% BOTTOM-RIGHT, EVEN CORE k LAYER #32, #80
%
if i==L && j==L && k~=2*L && mod(k,2)==0
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 +1 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2;
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -2*L^2 +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -2*L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2;
    %
end
    
%
% TOP-LEFT, EVEN CORE k LAYER
%
if i<L && j<L && k~=1 && k~=2*L && mod(k,2)==0
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 +1;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2 +L +1;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2;
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 +1;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2 +L +1;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2;
    %
end

%
% TOP-LEFT CORNER, ODD CORE k LAYER
%
if i==1 && j==1 && k~=1 && mod(k,2)~=0
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L +L^2 ;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1 +L;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -L +L^2 +L -1;
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L +L^2;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 -1 +L;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2 -L +L^2 +L -1;
    %
end

%
% BOTTOM-RIGHT, ODD k LAYER
%
if i>1 && j>1 && k~=1 && k~=2*L && mod(k,2)~=0
    %
    NN1_zpos1 = j+(i-1)*L+(k-1)*L*L +L^2 -L;
    NN1_zpos2 = j+(i-1)*L+(k-1)*L*L +L^2;
    NN1_zpos3 = j+(i-1)*L+(k-1)*L*L +L^2 -1;
    NN1_zpos4 = j+(i-1)*L+(k-1)*L*L +L^2 -L -1;
    %
    NN1_zneg1 = j+(i-1)*L+(k-1)*L*L -L^2 -L;
    NN1_zneg2 = j+(i-1)*L+(k-1)*L*L -L^2;
    NN1_zneg3 = j+(i-1)*L+(k-1)*L*L -L^2 -1;
    NN1_zneg4 = j+(i-1)*L+(k-1)*L*L -L^2 -L -1;
    %
end



%
% 2NN (SC - like)
%
% CORNERS [8]
%
if i == 1 && j == 1 && k == 1
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = L; %
    NN2_ypos = L^2-L+1; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = j+(i-1)*L+(k-1)*L*L + 2*L^2;
    NN2_zneg = 2*L^3-2*L^2+1; %
end
%
if i == 1 && j == L && k == 1
    NN2_xpos = 1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = L^2; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = 2*L^3-2*L^2+L; %
end
%
if i == L && j == 1 && k == 1
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = L^2; %
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = 1; %
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = 2*L^3-L^2-L+1; %
end
%
if i == L && j == L && k == 1
    NN2_xpos = L^2-L+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = L; %
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = 2*L^3-L^2; %
end
%
if i == 1 && j == 1 && k == 2*L % #113 
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = 2*L^3-L^2+L; %
    NN2_ypos = 2*L^3-L+1; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = L^2 + 1;
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
if i == 1 && j == L && k == 2*L % #116
    NN2_xpos = 2*L^3-L^2+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = 2*L^3; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = L^2 + L; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
if i == L && j == 1 && k == 2*L
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = 2*L^3; %
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = 2*L^3-L^2+1; %
    NN2_zpos = 2*L^2-L+1; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
if i == L && j == L && k == 2*L
    NN2_xpos = 2*L^3-L+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = 2*L^3-L^2+L; %
    NN2_zpos = 2*L^2; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
% EDGES [12]
%
if i == 1 && j~=1 && j~=L && k == 1 %#2
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = L^2-L+j; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = 2*L^3-2*L^2+j; %
end
%
if i == L && j~=1 && j~=L && k == 1 %#8
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j; %
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = j+(i-1)*L+(k-1)*L*L +2*L^3 -2*L^2; %
end
%
if i~=1 && i~=L && j == 1 && k == 1 %#4
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+L-1; %
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = 2*L^3-2*L^2+(i-1)*L+j;
end
%
if i~=1 && i~=L && j == L && k == 1 %#6
    NN2_xpos = j+(i-1)*L-L+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = 2*L^3-2*L^2+(i-1)*L+j;
end
%
if i == 1 && j == 1 && k~=1 && k~=2*L %#10
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = (k-1)*L^2+L; %
    NN2_ypos = (k-1)*L^2+L^2-L+1; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    %
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L +2*L^3 -2*L^2; 
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3; 
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i == 1 && j == L && k~=1 && k~=2*L %#12
    NN2_xpos = (k-1)*L^2+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = (k-1)*L^2+L^2; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i == L && j == 1 && k~=1 && k~=2*L %#16
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = (k-1)*L^2+L^2; %
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = (k-1)*L^2+1; %
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i == L && j == L && k~=1 && k~=2*L %#18
    NN2_xpos = (k-1)*L^2+L^2-L+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = (k-1)*L^2+L; %
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i == 1 && j~=1 && j~=L && k == 2*L % #114
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L +L^2 -L; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = L^2 + j; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
if i == L && j~=1 && j~=L && k == 2*L % #126 
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L -L^2 +L; %
    NN2_zpos = 2*L^2-L+j; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
if i~=1 && i~=L && j == 1 && k == 2*L  % #117
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L +L -1; %
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = L^2 + (i-1)*L+1; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
%
if i~=1 && i~=L && j == L && k == 2*L % #120
    NN2_xpos = j+(i-1)*L+(k-1)*L*L -L +1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = L^2 + i*L; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
% FACES [6]
%
if i~=1 && i~=L && j~=1 && j~=L && k == 1 %#5
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
    NN2_zneg = (L-1)*2*L^2+(i-1)*L+j; %
end
%
if i == 1 && j~=1 && j~=L && k~=1 && k~=2*L %#11
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = k*L^2-L+j; %
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i == L && j~=1 && j~=L && k~=1 && k~=2*L %#17
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = (k-1)*L^2+j; %
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i~=1 && i~=L && j == 1 && k~=1 && k~=2*L %#13
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    %NN2_xneg = (k-1)*L^2+L+i*L; % FIX = 15
    NN2_xneg = (k-1)*L^2+(i)*L; % TESTING
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i~=1 && i~=L && j == L && k~=1 && k~=2*L %#15
    NN2_xpos = (k-1)*L^2+(i-1)*L+1; %
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
%
if i~=1 && i~=L && j~=1 && j~=L && k == 2*L %#118
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    NN2_zpos = L^2 + (i-1)*L+j; %
    NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
end
%
% CORE
%
if i~=1 && i~=L && j~=1 && j~=L && k~=1 && k~=2*L
    NN2_xpos = j+(i-1)*L+(k-1)*L*L +1;
    NN2_xneg = j+(i-1)*L+(k-1)*L*L -1;
    NN2_ypos = j+(i-1)*L+(k-1)*L*L -L;
    NN2_yneg = j+(i-1)*L+(k-1)*L*L +L;
    if k ~= 2 && k ~= 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
    if k == 2
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2 +2*L^3;
    end
    %
    if k == 2*L-1
        NN2_zpos = j+(i-1)*L+(k-1)*L*L +2*L^2 -2*L^3;
        NN2_zneg = j+(i-1)*L+(k-1)*L*L -2*L^2;
    end
    %
end
