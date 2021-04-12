clear
close

% Validation for MPI new rps
% João Inácio, 29 Dec., 2020

L = 4;
N_SPINS = L * L;
skip = N_SPINS;

REP_hist = 10000;

REP_vals = [1000, 10000, 100000, 1000000];
idx_slice = 5;

run_time_avg = zeros(length(REP_vals), 1);

checkerboard = cell(length(REP_vals), 1);   % checkerboard - 2
slice = cell(length(REP_vals), 1);          % slice - 2*L
zerozero = cell(length(REP_vals), 1);       % zerozero - E, M = 0

mu = zeros(length(REP_vals), 3);
sigma = zeros(length(REP_vals), 3);
  
load("Data/JDOS_2D_" + L + "L_skip" + skip + "_nRPS_CPP.mat")

for idx_REP = 1:length(REP_vals)
    REP = REP_vals(idx_REP);
    
    for run = 1:1000
        checkerboard{idx_REP}(run) = JDOS{run, idx_REP}(N_SPINS + 1, N_SPINS / 2 + 1);
        slice{idx_REP}(run) = JDOS{run, idx_REP}(idx_slice, N_SPINS / 2 + 1);
        zerozero{idx_REP}(run) = JDOS{run, idx_REP}(N_SPINS / 2 + 1, N_SPINS / 2 + 1);
    end
    
    run_time_avg(idx_REP, 1) = mean(run_time(:, idx_REP));
    
    pd_chkbrd = fitdist(checkerboard{idx_REP}', 'normal');
    pd_slice = fitdist(slice{idx_REP}', 'normal');
    pd_zerozero = fitdist(zerozero{idx_REP}', 'normal');
    
    mu(idx_REP, 1) = pd_chkbrd.mu;
    mu(idx_REP, 2) = pd_slice.mu;
    mu(idx_REP, 3) = pd_zerozero.mu;
    
    sigma(idx_REP, 1) = pd_chkbrd.sigma;
    sigma(idx_REP, 2) = pd_slice.sigma;
    sigma(idx_REP, 3) = pd_zerozero.sigma;
end


if L == 4
    load("JDOS_exact_L4.mat")
    
    mean_error = zeros(1, length(REP_vals));
    
    load("Data/JDOS_2D_" + L + "L_skip" + skip + "_nRPS_CPP.mat")
    
    JDOS_error = cell(1000, length(REP_vals));
    error = zeros(1000, length(REP_vals));
    
    for idx_REP = 1:length(REP_vals)
        REP = REP_vals(idx_REP);
        
        for run = 1:1000
            JDOS_error{run, idx_REP} = abs(JDOS{run, idx_REP}(:, 1:N_SPINS / 2 + 1) - JDOS_exact(:, 1:N_SPINS / 2 + 1)) ./ JDOS_exact(:, 1:N_SPINS / 2 + 1);
            error(run, idx_REP) = sum(sum(JDOS_error{run, idx_REP}(~isnan(JDOS_error{run, idx_REP}))));
        end
        
        mean_error(1, idx_REP) = mean(error(:, idx_REP));
    end
        
    
    figure(1)
    subplot(1, 2, 1)
    plot(log10(REP_vals), mean_error(1, :), '-o')
    legend("Error - MPI", 'location', "best")
    title("Error vs log10(REP)")
    xlabel("log(REP)")
    
    subplot(1, 2, 2)
    plot(log10(REP_vals), run_time_avg(:, 1), '-o')
    legend("Run Time - MPI", 'location', "NW")
    title("Run Time vs long10(REP)")
    xlabel("log(REP)")

end

% figure(2)
% subplot(1, 2, 1)
% plot(REP_vals, run_time_avg(:, 1), '-o')
% legend("Run Time - MPI", 'location', "NW")
% title("Run Time vs REP")
% xlabel("REP")
% 
% subplot(1, 2, 2)
% plot(log10(REP_vals), run_time_avg(:, 1), '-o')
% legend("Run Time - MPI", 'location', "NW")
% title("Run Time vs long10(REP)")
% xlabel("log(REP)")

figure(3)
subplot(1, 3, 1)
% [plot_checkerboard, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 1), log10(REP_vals), mu{2}(:, 1), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_checkerboard(1).YLim = [0 1.5];
% plot_checkerboard(2).YLim = [0 400];
plot(log10(REP_vals), mu(:, 1), '-o')
legend("mu MPI", 'location', "best")
title("mu vs REP - checkerboard")
xlabel("log(REP)")

subplot(1, 3, 2)
% [plot_slice, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 2), log10(REP_vals), mu{2}(:, 2), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_slice(1).YLim = [0 1.5];
% plot_slice(2).YLim = [0 400];
plot(log10(REP_vals), mu(:, 2), '-o')
legend("mu MPI", 'location', "best")
title("mu vs log10(REP) - slice")
xlabel("log(REP)")

subplot(1, 3, 3)
% [plot_zerozero, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals),mu{1}(:, 3), log10(REP_vals), mu{2}(:, 3), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_zerozero(1).YLim = [0 1.5];
% plot_zerozero(2).YLim = [0 400];
plot(log10(REP_vals), mu(:, 3), '-o')
legend("mu MPI", 'location', "best")
title("mu vs log10(REP) - zerozero")
xlabel("log(REP)")

figure(4)
subplot(1, 3, 1)
plot(log10(REP_vals), sigma(:, 1), '-o')
legend("sigma MPI", 'location', "best")
title("sigma vs log10(REP) - checkerboard")
xlabel("log(REP)")

subplot(1, 3, 2)
plot(log10(REP_vals), sigma(:, 2), '-o')
legend("sigma MPI", 'location', "best")
title("sigma vs log10(REP) - slice")
xlabel("log(REP)")

subplot(1, 3, 3)
plot(log10(REP_vals), sigma(:, 3), '-o')
legend("sigma MPI", 'location', "best")
title("sigma vs log10(REP) - zerozero")
xlabel("log(REP)")

figure(5)
subplot(1,3,1)
h_chkbrd = histfit(checkerboard{REP_vals == REP_hist}, [], 'normal');
%     dataObjs_X = findobj(h_chkbrd,'-property','XData');
%     x1_chkbrd(:,idx_REP) = dataObjs_X(1).XData;
%     dataObjs_Y = findobj(h_chkbrd,'-property','YData');
%     y1_chkbrd(:,idx_REP) = dataObjs_Y(1).YData;
title("Checkerbord configuration - REP = " + REP_hist)

subplot(1,3,2)
h_slice = histfit(slice{REP_vals == REP_hist}, [], 'normal');
%     dataObjs_X = findobj(h_slice,'-property','XData');
%     x1_slice(:,idx_REP) = dataObjs_X(1).XData;
%     dataObjs_Y = findobj(h_slice,'-property','YData');
%     y1_slice(:,idx_REP) = dataObjs_Y(1).YData;
title("Slice configuration - REP = " + REP_hist)

subplot(1,3,3)
h_zerozero = histfit(zerozero{REP_vals == REP_hist}, [], 'normal');
%     dataObjs_X = findobj(h_zerozero,'-property','XData');
%     x1_zerozero(:,idx_REP) = dataObjs_X(1).XData;
%     dataObjs_Y = findobj(h_zerozero,'-property','YData');
%     y1_zerozero(:,idx_REP) = dataObjs_Y(1).YData;
title("E=0,M=0 configuration - REP = " + REP_hist)
















