function [nnxpos,nnxneg,nnypos,nnyneg,nnzpos,nnzneg]=function_NN_list_3D_SC(L,i,j,k)
%
% LATTICE CONVENTION:
% i=rows , j=columns , k=layers
% top row = lowest i value
% left column = lowest j value
% bottom layer = lowest k value
%
% CORNERS [8]
%
if i==1 && j==1 && k==1
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=L; %
    nnypos=L^2-L+1; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L^2+1; %
end
%
if i==1 && j==L && k==1
    nnxpos=1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=L^2; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L^2+L; %
end
%
if i==L && j==1 && k==1
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=L^2; %
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=1; %
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L+1; %
end
%
if i==L && j==L && k==1
    nnxpos=L^2-L+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=L; %
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3; %
end
%
if i==1 && j==1 && k==L
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=L^3-L^2+L; %
    nnypos=L^3-L+1; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=1; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==1 && j==L && k==L
    nnxpos=L^3-L^2+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=L^3; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=L; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==L && j==1 && k==L
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=L^3; %
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=L^3-L^2+1; %
    nnzpos=L^2-L+1; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==L && j==L && k==L
    nnxpos=L^3-L+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=L^3-L^2+L; %
    nnzpos=L^2; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
% EDGES [12]
%
if i==1 && j~=1 && j~=L && k==1 %#2
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=L^2-L+j; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L^2+j; %
end
%
if i==L && j~=1 && j~=L && k==1 %#8
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j; %
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L+j; %
end
%
if i~=1 && i~=L && j==1 && k==1 %#4
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+L-1; %
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L^2+(i-1)*L+j;
end
%
if i~=1 && i~=L && j==L && k==1 %#6
    nnxpos=j+(i-1)*L-L+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=L^3-L^2+(i-1)*L+j;
end
%
if i==1 && j==1 && k~=1 && k~=L %#10
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=(k-1)*L^2+L; %
    nnypos=(k-1)*L^2+L^2-L+1; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==1 && j==L && k~=1 && k~=L %#12
    nnxpos=(k-1)*L^2+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=(k-1)*L^2+L^2; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==L && j==1 && k~=1 && k~=L %#16
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=(k-1)*L^2+L^2; %
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=(k-1)*L^2+1; %
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==L && j==L && k~=1 && k~=L %#18
    nnxpos=(k-1)*L^2+L^2-L+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=(k-1)*L^2+L; %
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==1 && j~=1 && j~=L && k==L %#20
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=L^3-L+j; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==L && j~=1 && j~=L && k==L %#26
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=L^3-L^2+j; %
    nnzpos=L^2-L+j; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i~=1 && i~=L && j==1 && k==L %#22
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=L^3-L^2+i*L; %
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=(i-1)*L+1; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
%
if i~=1 && i~=L && j==L && k==L %#24
    nnxpos=L^3-L^2+(i-1)*L+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=i*L; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
% FACES [6]
%
if i~=1 && i~=L && j~=1 && j~=L && k==1 %#5
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=(L-1)*L^2+(i-1)*L+j; %
end
%
if i==1 && j~=1 && j~=L && k~=1 && k~=L %#11
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=k*L^2-L+j; %
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i==L && j~=1 && j~=L && k~=1 && k~=L %#17
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=(k-1)*L^2+j; %
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i~=1 && i~=L && j==1 && k~=1 && k~=L %#13
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    %nnxneg=(k-1)*L^2+L+i*L; % FIX = 15
    nnxneg=(k-1)*L^2+(i)*L; % TESTING
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i~=1 && i~=L && j==L && k~=1 && k~=L %#15
    nnxpos=(k-1)*L^2+(i-1)*L+1; %
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
if i~=1 && i~=L && j~=1 && j~=L && k==L %#23
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=(i-1)*L+j; %
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end
%
% CORE
%
if i~=1 && i~=L && j~=1 && j~=L && k~=1 && k~=L
    nnxpos=j+(i-1)*L+(k-1)*L*L +1;
    nnxneg=j+(i-1)*L+(k-1)*L*L -1;
    nnypos=j+(i-1)*L+(k-1)*L*L -L;
    nnyneg=j+(i-1)*L+(k-1)*L*L +L;
    nnzpos=j+(i-1)*L+(k-1)*L*L +L^2;
    nnzneg=j+(i-1)*L+(k-1)*L*L -L^2;
end

