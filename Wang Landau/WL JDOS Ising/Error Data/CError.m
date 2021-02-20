clear
close

% Computes the error for various JDOS with different f values
% João Inácio, Sep. 22, 2020

L = 4;
NSpins = L * L;

fVals = 1 + [1, 1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7,...
    1E-8, 1E-9, 1E-10, 1E-11, 1E-12, 1E-13, 1E-14, 1E-15];


load("JDOS_approx_L4.mat", "JDOS_approx_par")
load("JDOS_exact_L4.mat", "JDOS_exact")

% remove the runtime from the JDOS matrix, compute the averages for each
% value o f, and also computes the JDOS_error for each JDOS

runTime = zeros(1, length(fVals));
JDOS_approx = cell(1, length(fVals));

JDOS_error = cell(100, length(fVals));


for i = 1:length(fVals)
    JDOS_approx{i} = zeros(NSpins + 1);
    for it = 1:100
        JDOS_error{it, i} = abs(JDOS_approx_par{it, i}(:, 1:NSpins + 1) - JDOS_exact) ./ JDOS_exact;
        JDOS_approx{i} = JDOS_approx{i} + JDOS_approx_par{it, i}(:, 1:NSpins + 1);
        runTime(i) = runTime(i) + JDOS_approx_par{it, i}(1, NSpins + 2);
    end
    runTime(i) = runTime(i) / 100;
    JDOS_approx{i} = JDOS_approx{i} / 100;
end

% compute the error for each value of f

error = zeros(100, length(fVals));

for i = 1:length(fVals)
    for it = 1:100
        error(it, i) = sum(sum(JDOS_error{it, i}(~isnan(JDOS_error{it, i}))));
    end
end

mean_error = mean(error);
deviation = abs(mean_error - error);
err = max(deviation);

figure(1)
plotyy(abs(log10(1-fVals(4:end))), mean_error(4:end), abs(log10(1-fVals(4:end))), runTime(4:end))
legend("Error", "Runtime", 'location', "best")
title("Error and Runtime vs Modification Factor")
xlabel("f values")

figure(2)
errorbar(abs(log10(1-fVals(5:end))), mean_error(5:end), err(5:end))
title("Error vs Modification Factor with error bars")
ylabel("error")
xlabel("f values")

% get the f_idx th 100s JDOS to check if the error follows a normal
% distribution

% f_idx = 9;
% JDOS_error_10 = cell(1, 100);
% error_10 = zeros(1, 100);
% 
% for it = 1:100
%     JDOS_error_10{it} = abs(JDOS_approx_par{it, f_idx}(:, 1:NSpins + 1) - JDOS_exact) ./ JDOS_exact;
%     error_10(it) = sum(sum(JDOS_error_10{it}(~isnan(JDOS_error_10{it}))));
% end
% 
% figure(3)
% histfit(error_10, [], "normal")
% title("Distibution for the " + f_idx + " th value of f of the JDOS error")










