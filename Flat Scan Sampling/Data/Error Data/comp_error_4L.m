clear
close

% Comparison between the error of the new RPS with different REP values
% João Inácio, 18 Oct., 2020


load("JDOS_exact_L4.mat", "JDOS_exact");
runtime = load("runtime_2D_4L_nRPS_CPP.txt");

REP = [1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7];

JDOS_approx = cell(1, length(REP));
error = zeros(1, length(REP));

for i = 1:length(REP)
    JDOS = load("JDOS_2D_4L_10E" + i + "_nRPS_CPP.txt");
    
    JDOS(:, floor(length(JDOS) / 2) + 2:end) = JDOS(:, floor(length(JDOS) / 2):-1:1);
    
    JDOS_error = abs(JDOS - JDOS_exact) ./ JDOS_exact;
    error(i) =  sum(sum(JDOS_error(~isnan(JDOS_error))));
    
    JDOS_approx{1, i} = JDOS;
end

mean_error = mean(error);
deviation = abs(mean_error - error);
err = max(deviation);

figure(1)
plotyy(log10(REP), error, log10(REP), runtime)
legend("Error", "Runtime", 'location', "NE")
title("Error and Runtime vs log of REP value")
xlabel("log(REP)")

figure(2)
errorbar(log10(REP), error, deviation)
title("Error vs REP values with error bars")
ylabel("error")
xlabel("log(REP)")
