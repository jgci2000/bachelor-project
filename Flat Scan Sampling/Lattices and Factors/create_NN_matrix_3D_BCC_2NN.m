clear
close all
%
L = 16;
%
NN_1 = 8;
NN_2 = 6;
%
NN_matrix = nan(2*L^3, NN_1 + NN_2);
%
for i = 1:L
    for j = 1:L
        for k = 1:(2*L)
            %
            [ NN_matrix(j+(i-1)*L+(k-1)*L*L, 1), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 2), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 3), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 4), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 5), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 6), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 7), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 8), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 9), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 10), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 11), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 12), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 13), ...
                NN_matrix(j+(i-1)*L+(k-1)*L*L, 14), ...
                ] = function_NN_2NN_list_BCC(L,i,j,k);
            %
        end
    end
end
%

%
% CHECK OUT Of BOUNDS, REPETITIONS AND SELF NEIGHBOURS (FOR BCC bonds only 1:8)
%
disp('checking BCC (1:8) NN list for values out of bounds')
%
for n = 1:(2*L^3)
    %
    if any(NN_matrix(n, 1:8) < 1) || any(NN_matrix(n, 1:8) > 2*L^3)
        %
        disp(['atom number ', int2str(n), ' has BCC NN values out of bounds! FIX!'])
        %
    end
    %
end

disp('checking BCC (1:8) NN list for repetitions')
%
for n = 1:(2*L^3)
    %
    if length(unique(NN_matrix(n, 1:8))) ~= 8
        %
        disp(['atom number ', int2str(n), ' has BCC repeated neighbours! FIX!'])
        %
    end
    %
end
%
for N_direction = 1:8
    %
    if length(unique(NN_matrix(:, N_direction))) ~= 2*L^3
        %
        disp(['neighbour direction ', int2str(N_direction), ' has BCC repeated neighbours! FIX!'])
        %
    end
    %
end
%
disp('checking BCC (1:8) NN list for self-neighbours')
%
for n = 1:(2*L^3)
    %
    if any(NN_matrix(n, 1:8) == n)
        %
        disp(['atom number ', int2str(n), ' has a BCC self-neighbour! FIX!'])
        %
    end
    %
end
%
% SC
%
disp('checking SC (9:14) NN list for values out of bounds')
%
for n = 1:(2*L^3)
    %
    if any(NN_matrix(n, 9:14) < 1) || any(NN_matrix(n, 9:14) > 2*L^3)
        %
        disp(['atom number ', int2str(n), ' has SC NN values out of bounds! FIX!'])
        %
    end
    %
end

disp('checking SC (9:14) NN list for repetitions')
%
for n = 1:(2*L^3)
    %
    if length(unique(NN_matrix(n, 9:14))) ~= 6
        %
        disp(['atom number ', int2str(n), ' has SC repeated neighbours! FIX!'])
        %
    end
    %
end
%
for N_direction = 9:14
    %
    if length(unique(NN_matrix(:, N_direction))) ~= 2*L^3
        %
        disp(['neighbour direction ', int2str(N_direction), ' has SC repeated neighbours! FIX!'])
        %
    end
    %
end
%
disp('checking SC (1:8) NN list for self-neighbours')
%
for n = 1:(2*L^3)
    %
    if any(NN_matrix(n, 9:14) == n)
        %
        disp(['atom number ', int2str(n), ' has a SC self-neighbour! FIX!'])
        %
    end
    %
end
%
% BOTH BCC AND SC CHECKS
%
disp('checking both BCC and SC (1:14) NN list for repetitions')
%
for n = 1:(2*L^3)
    %
    if length(unique(NN_matrix(n, 1:14))) ~= 14
        %
        disp(['atom number ', int2str(n), ' has repeated neighbours between BCC and SC! FIX!'])
        %
    end
    %
end
%
% CHECK NEIGHBOUR RECIPROCIPITY
%
disp('checking both BCC and SC (1:14) NN list for reciprocipity')
%
for n = 1:(2*L^3)
    %
    if n ~= NN_matrix(NN_matrix(n, 1), 7)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 1)), ' on directions 1 and 7! FIX!'])
        %
    end
    %
    if n ~= NN_matrix(NN_matrix(n, 2), 8)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 2)), ' on directions 2 and 8! FIX!'])
        %
    end
    %
    if n ~= NN_matrix(NN_matrix(n, 3), 5)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 3)), ' on directions 3 and 5! FIX!'])
        %
    end
    %
    if n ~= NN_matrix(NN_matrix(n, 4), 6)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 1)), ' on directions 4 and 6! FIX!'])
        %
    end
    %
    if n ~= NN_matrix(NN_matrix(n, 9), 10)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 1)), ' on directions 9 and 10! FIX!'])
        %
    end
    %
    if n ~= NN_matrix(NN_matrix(n, 11), 12)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 1)), ' on directions 11 and 12! FIX!'])
        %
    end
    %
    if n ~= NN_matrix(NN_matrix(n, 13), 14)
        %
        disp(['no reciprocity between atom ', int2str(n), ' and atom ', int2str(NN_matrix(n, 1)), ' on directions 13 and 14! FIX!'])
        %
    end
end
%
NN_table = NN_matrix;
eval(['save neighbour_table_3D_BCC_2NN_L', int2str(L), ' NN_table -v7.3'])