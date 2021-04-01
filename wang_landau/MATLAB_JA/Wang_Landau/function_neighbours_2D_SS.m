function [nnxpos, nnxneg, nnypos, nnyneg] = function_neighbours_2D_SS(L)
%
% LATTICE AND NEIGHBOUR TABLES
%
lattice = repmat(1:L,L,1);
%
for index = 1:L
    %
    lattice(index,:) = lattice(index,:) + (index-1)*L;
    %
end
%
nnxpos = circshift(lattice,[0,-1]);
nnxneg = circshift(lattice,[0,1]);
nnypos = circshift(lattice,[1,0]);
nnyneg = circshift(lattice,[-1,0]);
%
nnxpos = reshape(nnxpos', [], 1);
nnxneg = reshape(nnxneg', [], 1);
nnypos = reshape(nnypos', [], 1);
nnyneg = reshape(nnyneg', [], 1);
