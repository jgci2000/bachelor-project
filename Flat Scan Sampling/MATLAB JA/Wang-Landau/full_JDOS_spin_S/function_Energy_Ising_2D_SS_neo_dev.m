function U = function_Energy_Ising_2D_SS_neo_dev(L, S_vector, NN_table)
%
% ENERGY CALCULATION
%
U=0;
%
for i = 1:L
    %
    for j = 1:L
        %
        for a = 1:length(NN_table(1,:))
            %
            U = U - 1/2*S_vector(j+(i-1)*L,1)*(S_vector(NN_table(j+(i-1)*L,a)));
            %
        end
        %
    end
    %
end
%
end