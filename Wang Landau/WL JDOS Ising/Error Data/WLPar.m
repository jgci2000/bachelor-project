clear
close

% Computes the JDOS from the WL sampling for various f values
% João Inácio, Sep. 18, 2020

L = 4;
fVals = 1 + [1, 1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7,...
    1E-8, 1E-9, 1E-10, 1E-11, 1E-12, 1E-13, 1E-14, 1E-15];


JDOS_approx_par = cell(100, length(fVals));

t_init = tic;
for NF = 1:length(fVals)
    fMin = fVals(NF);
    
    t_par = tic;
    parfor k = 1:100
        JDOS_approx_par{k, NF} = WLPar(L, fMin);
    end
    end_t = toc(t_par);
    fprintf("100 times took %.3f seconds to run with f = %.15f\n", end_t, fMin);
end
toc(t_init)


