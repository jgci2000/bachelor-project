function U=function_Energy_Ising_2D_SS(L,S_vector)
%
% VECTOR NEIGHBOUR TABLE
%
nnxpos=nan(L^2,1);
nnxneg=nan(L^2,1);
nnypos=nan(L^2,1);
nnyneg=nan(L^2,1);
%
for i=1:L
    for j=1:L
            [nnxpos(j+(i-1)*L),nnxneg(j+(i-1)*L),nnypos(j+(i-1)*L),nnyneg(j+(i-1)*L)]=function_NN_list_2D_SS(L,i,j);
    end
end
%
% ENERGY CALCULATION
%
U=0;
%
for i=1:L
    for j=1:L
        U = U - 1/2*S_vector(j+(i-1)*L,1)*(S_vector(nnxpos(j+(i-1)*L))+S_vector(nnxneg(j+(i-1)*L))+S_vector(nnypos(j+(i-1)*L))+S_vector(nnyneg(j+(i-1)*L)));
    end
end
%