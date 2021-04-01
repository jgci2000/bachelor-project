clear
close

% Compute the thermodynamics for 2D Ising Model from the JDOS
% João Inácio, Sep. 7, 2020

LVals = [2, 4, 8];

for k = 1:length(LVals)
L = LVals(k);

load("JDOS_2D_" + L + "L_WL.mat", "g");
gCPP = load("JDOS_2D_" + L + "L_WL_CPP.txt");

NSpins = L^2;     % Number of spins

% All of the possible energies
energies = 4 * (0:NSpins) - 2 * NSpins;
NE = length(energies);

% Magnetizations
magnetizations = - NSpins:2:NSpins;
NM = length(magnetizations);

% Get the DOS
gE = zeros(1, NE);
for i = 1:NE
    gE(1, i) = sum(g(i, :));
end

gECPP = zeros(1, L);
for i = 1:NE
    gECPP(1, i) = sum(gCPP(i, :));
end

% Temperature vector
% Temperatures below 0.75 give - Inf or NaN values on thermodynamic
% properties
T = 0.75:0.1:20;
NT = length(T);
kB = 1;
beta = zeros(1, NT);

Z = zeros(1, NT);
U = zeros(1, NT);
C = zeros(1, NT);
F = zeros(1, NT);
S = zeros(1, NT);
M = zeros(1, NT);

ZCPP = zeros(1, NT);
UCPP = zeros(1, NT);
CCPP = zeros(1, NT);
FCPP = zeros(1, NT);
SCPP = zeros(1, NT);
MCPP = zeros(1, NT);

for i = 1:NT
    beta(i) = - 1 / (T(i) * kB);
    
    % Partition Function
    Z(i) = sum(gE .* exp(energies * beta(i)));
    ZCPP(i) = sum(gECPP .* exp(energies * beta(i)));
    
    % Compute the Internal Energy, U(T)
    U(i) = sum(energies .* gE .* exp(energies .* beta(i))) / Z(i);
    UCPP(i) = sum(energies .* gECPP .* exp(energies * beta(i))) / ZCPP(i);
    
    % Compute the Specific Heat, C(T)
    avgE = U(i);
    avgE2 = sum(energies.^2 .* gE .* exp(energies .* beta(i))) / Z(i);
    C(i) = (avgE2 - avgE^2) / (kB * T(i)^2);
    avgECPP = UCPP(i);
    avgE2CPP = sum(energies.^2 .* gECPP .* exp(energies * beta(i))) / ZCPP(i);
    CCPP(i) = (avgE2CPP - avgECPP^2) / (kB * T(i)^2);
    
    % Compute the Free Energy, F(T)
    F(i) = -kB * T(i) * log(Z(i));
    FCPP(i) = -kB * T(i) * log(ZCPP(i));
    
    % Compute the Entropy, S(T)
    S(i) = (U(i) - F(i)) / T(i);
    SCPP(i) = (UCPP(i) - FCPP(i)) / T(i);

    % Compute the Mean Magnetization, <|M|>
    for idxE = 1:NE
        for idxM = 1:NM
            M(i) = M(i) + sum(g(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / Z(i);
            MCPP(i) = MCPP(i) + sum(gCPP(idxE, idxM) * abs(magnetizations(idxM)) * exp(energies(idxE) * beta(i))) / ZCPP(i);
        end
    end
end

% Approximations of Tc
TcExact = 2.26918531421;

Tc_C_Mat(k) = T(C == max(C));
% Tc_C_Cpp(k) = T(CCPP == max(CCPP));

figure(2)
plot(T, M/NSpins)
hold on
legend("L = 2", "L = 4", "L = 8");
xlabel("Temperature")
ylabel("Mean Magnetization")

figure(3)
subplot(2, 2, 1)
plot(T, U)
hold on
legend("L = 2", "L = 4", "L = 8", 'Location', "SE")
xlabel("Temperature")
ylabel("Internal Energy")

subplot(2, 2, 2)
plot(T, C)
hold on
legend("L = 2", "L = 4", "L = 8")
xlabel("Temperature")
ylabel("Specific Heat")

subplot(2, 2, 3)
plot(T, F)
hold on
legend("L = 2", "L = 4", "L = 8", 'Location', "SW")
xlabel("Temperature")
ylabel("Free Energy")

subplot(2, 2, 4)
plot(T, S)
hold on
legend("L = 2", "L = 4", "L = 8", 'Location', "NW")
xlabel("Temperature")
ylabel("Entropy")

sgtitle("Comparison of thermodynamic quantities from 3 Ising model lattices")
end


figure(4)
subplot(2, 3, 1)
plot(T, U)
hold on
plot(T, UCPP)
legend("Matlab", "C++", 'Location', "SE")
xlabel("Temperature")
ylabel("Internal Energy for L = 8")

subplot(2, 3, 2)
plot(T, C)
hold on
plot(T, CCPP)
legend("Matlab", "C++")
xlabel("Temperature")
ylabel("Specific Heat for L = 8")

subplot(2, 3, 3)
plot(T, F)
hold on
plot(T, FCPP)
legend("Matlab", "C++", 'Location', "SW")
xlabel("Temperature")
ylabel("Free Energy for L = 8")

subplot(2, 3, 4)
plot(T, S)
hold on
plot(T, SCPP)
legend("Matlab", "C++", 'Location', "NW")
xlabel("Temperature")
ylabel("Entropy for L = 8")

subplot(2, 3, 5)
plot(T, M/NSpins)
hold on
plot(T, MCPP/NSpins)
legend("MATLAB", "C++")
xlabel("Temperature")
ylabel("Mean Magnetization for L = 8")

subplot(2, 3, 6)


sgtitle("MATLAB vs C++ Ising model 2D 8 lattice")

% figure(5)
% Tc_C_Mat_Interp = interp1(Tc_C_Mat, [2, 4, 8], 2:32, 'linear');
% 
% plot([2, 4, 8], Tc_C_Mat, 'ko-')
% hold on
% % plot([2, 4, 8], Tc_C_Cpp, 'ro-')
% plot([2, 4, 8], [TcExact, TcExact, TcExact], 'b-')

