function [nnxpos,nnxneg,nnypos,nnyneg]=function_NN_list_2D_SS(L,i,j)
%
% LATTICE CONVENTION:
% i=rows , j=columns
% top row = lowest i value
% left column = lowest j value
%
% CORNERS [4]
%
if i==1 && j==1 % #1
    nnxpos=j+(i-1)*L +1;
    nnxneg=L; %
    nnypos=L^2-L +1 ; %
    nnyneg=j+(i-1)*L +L;
end
%
if i==1 && j==L % #3
    nnxpos=1; %
    nnxneg=j+(i-1)*L -1;
    nnypos=L^2; %
    nnyneg=j+(i-1)*L +L;
end
%
if i==L && j==1 % #7
    nnxpos=j+(i-1)*L +1;
    nnxneg=L^2; %
    nnypos=j+(i-1)*L -L;
    nnyneg=1; %
end
%
if i==L && j==L % #9
    nnxpos=L^2-L+1; %
    nnxneg=j+(i-1)*L -1;
    nnypos=j+(i-1)*L -L;
    nnyneg=L; %
end
%
% EDGES [4]
%
if i==1 && j~=1 && j~=L % #2
    nnxpos=j+(i-1)*L +1;
    nnxneg=j+(i-1)*L -1;
    nnypos=L^2-L+j; %
    nnyneg=j+(i-1)*L +L;
end
%
if j==1 && i~=1 && i~=L % #4
    nnxpos=j+(i-1)*L +1;
    nnxneg=i*L; %
    nnypos=j+(i-1)*L -L;
    nnyneg=j+(i-1)*L +L;
end
%
if j==L && i~=1 && i~=L % #6
    nnxpos=(i-1)*L +1; %
    nnxneg=j+(i-1)*L -1;
    nnypos=j+(i-1)*L -L;
    nnyneg=j+(i-1)*L +L;
end
%
if i==L && j~=1 && j~=L % #8
    nnxpos=j+(i-1)*L +1;
    nnxneg=j+(i-1)*L -1;
    nnypos=j+(i-1)*L -L;
    nnyneg=j; %
end
%
% CORE
%
if i~=1 && i~=L && j~=1 && j~=L
    nnxpos=j+(i-1)*L +1;
    nnxneg=j+(i-1)*L -1;
    nnypos=j+(i-1)*L -L;
    nnyneg=j+(i-1)*L +L;
end
