clear
close

load("JDOS_approx_L4.mat")
load("JDOS_exact_L4.mat")

JDOS_approx(:, 1:9) = JDOS_approx_par{1, 9}(:, 1:9);


sumExact = sum(sum(JDOS_exact(:, 1:9)));
sumApprox = sum(sum(JDOS_approx(:, 1:9)));

JDOS_error(:, 1:9) = abs(JDOS_approx(:, 1:9) - JDOS_exact(:, 1:9)) ./ JDOS_exact(:, 1:9);    

sum(sum(JDOS_error(~isnan(JDOS_error))))
