clear
close
clc

% Compute the thermodynamics for 2D Ising Model from the DOS
% João Inácio, Sep. 1, 2020

LVals = [2, 4, 8, 16];
Tc_C_Mat = zeros(1, length(LVals));
Tc_C_Cpp = zeros(1, length(LVals));

for k = 1:length(LVals)
L = LVals(k);
    
load("DOS_" + L + "L_WL.mat", "g");
gCPP = load("DOS_" + L + "L_WL_CPP.txt");

% L = 2;              % Number of spins for each side of the lattice.
NSpins = L * L;     % Number of spins.

% All of the possible energies. EMin + 4 and EMax -4 are not possible.
NE = NSpins - 1;
energies = zeros(1, NE);
temp = 4 * (0:NSpins) - 2 * NSpins;
energies(1) = temp(1);
energies(2:length(temp) - 3) = temp(3:end-2);
energies(length(temp) - 2) = temp(end);

% Temperature vector
% Temperatures below 0.75 give Inf or NaN values on thermodynamic
% properties
T = 0.75:0.01:8;
NT = length(T);
kB = 1;
beta = zeros(1, NT);
Z = zeros(1, NT);
U = zeros(1, NT);
C = zeros(1, NT);
F = zeros(1, NT);
S = zeros(1, NT);

ZCPP = zeros(1, NT);
UCPP = zeros(1, NT);
CCPP = zeros(1, NT);
FCPP = zeros(1, NT);
SCPP = zeros(1, NT);

for i = 1:NT
    beta(i) = - 1 / (T(i) * kB);
    
    % Partition Function
    Z(i) = sum(g .* exp(energies .* beta(i)));
    ZCPP(i) = sum(gCPP .* exp(energies .* beta(i)));
    
    % Compute the Internal Energy, U(T)
    U(i) = sum(energies .* g .* exp(energies .* beta(i))) ./ Z(i);
    UCPP(i) = sum(energies .* gCPP .* exp(energies .* beta(i))) ./ ZCPP(i);
    
    % Compute the Specific Heat, C(T)
    avgE = U(i);
    avgE2 = sum(energies.^2 .* g .* exp(energies .* beta(i))) ./ Z(i);
    C(i) = (avgE2 - avgE^2) / (kB * T(i)^2);
    avgECPP = UCPP(i);
    avgE2CPP = sum(energies.^2 .* gCPP .* exp(energies .* beta(i))) ./ ZCPP(i);
    CCPP(i) = (avgE2CPP - avgECPP^2) / (kB * T(i)^2);
    
    % Compute the Free Energy, F(T)
    F(i) = -kB .* T(i) .* log(Z(i));
    FCPP(i) = -kB .* T(i) .* log(ZCPP(i));
    
    % Compute the Entropy, S(T)
    S(i) = (U(i) - F(i)) / T(i);
    SCPP(i) = (UCPP(i) - FCPP(i)) / T(i);
end

% Approximations of Tc
TcExact = 2.26918531421;

Tc_C_Mat(k) = T(C == max(C));
Tc_C_Cpp(k) = T(CCPP == max(CCPP));

figure(3)
subplot(2, 2, 1)
plot(T, U)
hold on
legend("L = 2", "L = 4", "L = 8", "L = 16", 'Location', "SE")
xlabel("Temperature")
ylabel("Internal Energy")

subplot(2, 2, 2)
plot(T, C)
hold on
legend("L = 2", "L = 4", "L = 8", "L = 16")
xlabel("Temperature")
ylabel("Specific Heat")

subplot(2, 2, 3)
plot(T, F)
hold on
legend("L = 2", "L = 4", "L = 8", "L = 16", 'Location', "SW")
xlabel("Temperature")
ylabel("Free Energy")

subplot(2, 2, 4)
plot(T, S)
hold on
legend("L = 2", "L = 4", "L = 8", "L = 16", 'Location', "NW")
xlabel("Temperature")
ylabel("Entropy")

sgtitle("Comparison of thermodynamic quantities from 4 Ising model lattices")
end

figure(4)
subplot(2, 2, 1)
plot(T, U)
hold on
plot(T, UCPP)
legend("Matlab", "C++", 'Location', "SE")
xlabel("Temperature")
ylabel("Internal Energy for L = 16")

subplot(2, 2, 2)
plot(T, C)
hold on
plot(T, CCPP)
legend("Matlab", "C++")
xlabel("Temperature")
ylabel("Specific Heat for L = 16")

subplot(2, 2, 3)
plot(T, F)
hold on
plot(T, FCPP)
legend("Matlab", "C++", 'Location', "SW")
xlabel("Temperature")
ylabel("Free Energy for L = 16")

subplot(2, 2, 4)
plot(T, S)
hold on
plot(T, SCPP)
legend("Matlab", "C++", 'Location', "NW")
xlabel("Temperature")
ylabel("Entropy for L = 16")

sgtitle("MATLAB vs C++ Ising model 2D 16 lattice")

figure(5)

Tc_C_Mat_Interp = interp1(Tc_C_Mat, [2, 4, 8, 16], 2:32, 'linear');

plot([2, 4, 8, 16], Tc_C_Mat, 'ko-')
hold on
plot([2, 4, 8, 16], Tc_C_Cpp, 'ro-')
plot([2, 4, 8, 16], [TcExact, TcExact, TcExact, TcExact], 'b-')

