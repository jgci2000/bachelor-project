function U = function_Energy_Ising_3D_SC(L,S_vector)
%
% PRELIMINARY CALCULATION OF NEIGHBOUR TABLES
%
N_atm = L^3;
%
x_p=nan(N_atm,1);
x_n=nan(N_atm,1);
y_p=nan(N_atm,1);
y_n=nan(N_atm,1);
z_p=nan(N_atm,1);
z_n=nan(N_atm,1);
%
for i=1:L
    for j=1:L
        for k=1:L
            [x_p(j+(i-1)*L+(k-1)*L^2),x_n(j+(i-1)*L+(k-1)*L^2),y_p(j+(i-1)*L+(k-1)*L^2),y_n(j+(i-1)*L+(k-1)*L^2),z_p(j+(i-1)*L+(k-1)*L^2),z_n(j+(i-1)*L+(k-1)*L^2)]=function_NN_list_3D_SC(L,i,j,k);
        end
    end
end
%
% ENERGY CALCULATION
%
U=0;
%
for i=1:L
    %
    for j=1:L
        %
        for k=1:L
            %
            U = U - 1/2*(...
                S_vector(j+(i-1)*L+(k-1)*L^2,1) .* S_vector(x_p(j+(i-1)*L+(k-1)*L^2)) + ...
                S_vector(j+(i-1)*L+(k-1)*L^2,1) .* S_vector(x_n(j+(i-1)*L+(k-1)*L^2)) + ...
                S_vector(j+(i-1)*L+(k-1)*L^2,1) .* S_vector(y_p(j+(i-1)*L+(k-1)*L^2)) + ...
                S_vector(j+(i-1)*L+(k-1)*L^2,1) .* S_vector(y_n(j+(i-1)*L+(k-1)*L^2)) + ...
                S_vector(j+(i-1)*L+(k-1)*L^2,1) .* S_vector(z_p(j+(i-1)*L+(k-1)*L^2)) + ...
                S_vector(j+(i-1)*L+(k-1)*L^2,1) .* S_vector(z_n(j+(i-1)*L+(k-1)*L^2)) ) ;%-H*S_vector(j+(i-1)*L,1);
            %
        end
        %
    end
    %
end
%
%U=U/(L^2);