clear
close

% Computes the M vs T curve for the 3D Heisenberg model with the Metropolis
% algorithm with the help of parallelization
% João Inácio, Sep. 24, 2020

L = 16;
it_max = 1e6;

T_vals = 0:0.1:5;

M = cell(length(T_vals), 4);

t_init = tic;
for nT = 1:length(T_vals)
    T = T_vals(nT);
    
    t_par = tic;
    parfor k = 1:4
        M{nT, k} = Heisenberg(L, T, it_max);
    end
    t_end = toc(t_par);
    fprintf("Time elapsed: %.3fs for T = %.2f.\n", t_end, T)
end
toc(t_init)

