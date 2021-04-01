clear
close


% Plots M vs T from the parallellization data for the Heisenberg model
% João Inácio, Set. 24, 2020

load("M_16L_4x_010", "M")

L = 16;
N = L ^ 3;
it_max = 1e6;

T = 0:0.1:5;

mean_M = zeros(1, length(T));
mean_M_4 = zeros(length(T), 4);

for nT = 1:length(T)
    for i = 1:4
        mean_M_4(nT, i) = mean(abs(M{nT, i}((end - (it_max / 10)):end))) / N;
    end
    mean_M(nT) = mean(mean_M_4(nT, :));
end

deviation = abs(mean_M - mean_M_4.');
err = max(deviation);

figure(1)
errorbar(T, mean_M, err);
title("Magnetization vs T for " + L + " 3D lattice Heisenberg model with error bars")
xlabel("Temperature")
ylabel("Mean Magnetization per Spin")







