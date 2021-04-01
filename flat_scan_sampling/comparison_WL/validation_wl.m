clear
close

% Validation for MPI WL
% João Inácio, 29 Dec., 2020

L = 4;
N_SPINS = L * L;
flatness = 0.9;

exp_hist = 8;

exp_f_final_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
idx_slice = 5;

run_time_avg = zeros(length(exp_f_final_vals), 1);

checkerboard = cell(length(exp_f_final_vals), 1);   % checkerboard - 2
slice = cell(length(exp_f_final_vals), 1);          % slice - 2*L
zerozero = cell(length(exp_f_final_vals), 1);       % zerozero - E, M = 0

mu = zeros(length(exp_f_final_vals), 3);
sigma = zeros(length(exp_f_final_vals), 3);
  
load("Data/JDOS_2D_" + L + "L_flatness" + (flatness * 100) + "_WL_CPP.mat")

for idx_exp = 1:length(exp_f_final_vals)
    exp_f_final = exp_f_final_vals(idx_exp);
    
    for run = 1:1000
        checkerboard{idx_exp}(run) = JDOS{run, idx_exp}(N_SPINS + 1, N_SPINS / 2 + 1);
        slice{idx_exp}(run) = JDOS{run, idx_exp}(idx_slice, N_SPINS / 2 + 1);
        zerozero{idx_exp}(run) = JDOS{run, idx_exp}(N_SPINS / 2 + 1, N_SPINS / 2 + 1);
    end
    
    run_time_avg(idx_exp, 1) = mean(run_time(:, idx_exp));
    
    pd_chkbrd = fitdist(checkerboard{idx_exp}', 'normal');
    pd_slice = fitdist(slice{idx_exp}', 'normal');
    pd_zerozero = fitdist(zerozero{idx_exp}', 'normal');
    
    mu(idx_exp, 1) = pd_chkbrd.mu;
    mu(idx_exp, 2) = pd_slice.mu;
    mu(idx_exp, 3) = pd_zerozero.mu;
    
    sigma(idx_exp, 1) = pd_chkbrd.sigma;
    sigma(idx_exp, 2) = pd_slice.sigma;
    sigma(idx_exp, 3) = pd_zerozero.sigma;
end


if L == 4
    load("JDOS_exact_L4.mat")
    
    mean_error = zeros(1, length(exp_f_final_vals));
    
    load("Data/JDOS_2D_" + L + "L_flatness" + (flatness * 100) + "_WL_CPP.mat")
    
    JDOS_error = cell(1000, length(exp_f_final_vals));
    error = zeros(1000, length(exp_f_final_vals));
    
    for idx_exp = 1:length(exp_f_final_vals)
        exp_f_final = exp_f_final_vals(idx_exp);
        
        for run = 1:1000
            JDOS_error{run, idx_exp} = abs(JDOS{run, idx_exp}(:, 1:N_SPINS / 2 + 1) - JDOS_exact(:, 1:N_SPINS / 2 + 1)) ./ JDOS_exact(:, 1:N_SPINS / 2 + 1);
            error(run, idx_exp) = sum(sum(JDOS_error{run, idx_exp}(~isnan(JDOS_error{run, idx_exp}))));
        end
        
        mean_error(1, idx_exp) = mean(error(:, idx_exp));
    end
        
    
    figure(1)
    subplot(1, 2, 1)
    plot(exp_f_final_vals(3:end), mean_error(1, 3:end), '-o')
    legend("Error", 'location', "best")
    title("Error vs log10(f final)")
    xlabel("log(f final)")
    ylim([0 8])
    
    subplot(1, 2, 2)
    plot(exp_f_final_vals(3:end), run_time_avg(3:end, 1), '-o')
    legend("Run Time", 'location', "NW")
    title("Run Time vs long10(f final)")
    xlabel("log(f final)")
end

% figure(2)
% subplot(1, 2, 1)
% plot(1 + exp(- exp_f_final_vals(3:end)), run_time_avg(3:end, 1), '-o')
% legend("Run Time", 'location', "NW")
% title("Run Time vs f final")
% xlabel("f final")
% 
% subplot(1, 2, 2)
% plot(exp_f_final_vals(3:end), run_time_avg(3:end, 1), '-o')
% legend("Run Time", 'location', "NW")
% title("Run Time vs long10(f final)")
% xlabel("log(f final)")

figure(3)
subplot(1, 3, 1)
% [plot_checkerboard, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 1), log10(REP_vals), mu{2}(:, 1), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_checkerboard(1).YLim = [0 1.5];
% plot_checkerboard(2).YLim = [0 400];
plot(exp_f_final_vals(3:end), mu(3:end, 1), '-o')
legend("mu", 'location', "best")
title("mu vs log10(f final) - checkerboard")
xlabel("log(f final)")

subplot(1, 3, 2)
% [plot_slice, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 2), log10(REP_vals), mu{2}(:, 2), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_slice(1).YLim = [0 1.5];
% plot_slice(2).YLim = [0 400];
plot(exp_f_final_vals(3:end), mu(3:end, 2), '-o')
legend("mu", 'location', "best")
title("mu vs log10(f final) - slice")
xlabel("log(f final)")

subplot(1, 3, 3)
% [plot_zerozero, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals),mu{1}(:, 3), log10(REP_vals), mu{2}(:, 3), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_zerozero(1).YLim = [0 1.5];
% plot_zerozero(2).YLim = [0 400];
plot(exp_f_final_vals(3:end), mu(3:end, 3), '-o')
legend("mu", 'location', "best")
title("mu vs log10(f final) - zerozero")
xlabel("log(f final)")

figure(4)
subplot(1, 3, 1)
plot(exp_f_final_vals(3:end), sigma(3:end, 1), '-o')
legend("sigma", 'location', "best")
title("sigma vs log10(f final) - checkerboard")
xlabel("log(f final)")

subplot(1, 3, 2)
plot(exp_f_final_vals(3:end), sigma(3:end, 2), '-o')
legend("sigma", 'location', "best")
title("sigma vs log10(f final) - slice")
xlabel("log(f final)")

subplot(1, 3, 3)
plot(exp_f_final_vals(3:end), sigma(3:end, 3), '-o')
legend("sigma", 'location', "best")
title("sigma vs log10(f final) - zerozero")
xlabel("log(f final)")

figure(5)
subplot(1,3,1)
h_chkbrd = histfit(checkerboard{exp_f_final_vals == exp_hist}, [], 'normal');
%     dataObjs_X = findobj(h_chkbrd,'-property','XData');
%     x1_chkbrd(:,idx_REP) = dataObjs_X(1).XData;
%     dataObjs_Y = findobj(h_chkbrd,'-property','YData');
%     y1_chkbrd(:,idx_REP) = dataObjs_Y(1).YData;
title("Checkerbord configuration - log(f final) = " + exp_hist)

subplot(1,3,2)
h_slice = histfit(slice{exp_f_final_vals == exp_hist}, [], 'normal');
%     dataObjs_X = findobj(h_slice,'-property','XData');
%     x1_slice(:,idx_REP) = dataObjs_X(1).XData;
%     dataObjs_Y = findobj(h_slice,'-property','YData');
%     y1_slice(:,idx_REP) = dataObjs_Y(1).YData;
title("Slice configuration - log(f final) = " + exp_hist)

subplot(1,3,3)
h_zerozero = histfit(zerozero{exp_f_final_vals == exp_hist}, [], 'normal');
%     dataObjs_X = findobj(h_zerozero,'-property','XData');
%     x1_zerozero(:,idx_REP) = dataObjs_X(1).XData;
%     dataObjs_Y = findobj(h_zerozero,'-property','YData');
%     y1_zerozero(:,idx_REP) = dataObjs_Y(1).YData;
title("E=0,M=0 configuration - log(f final) = " + exp_hist)
















