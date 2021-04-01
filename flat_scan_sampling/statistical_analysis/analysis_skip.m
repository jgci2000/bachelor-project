clear
close

% Error analysis for the CPP results
% João Inácio, 12 Nov., 2020

% Change REP_hist to change which REP value is ploted in the histogram
% plots

L = 8;
N_SPINS = L * L;
type_vals = ["Par", "Serial"]; % Serial/Par

skip_hist = 64;

if L == 4
    skip_vals = [0 2 4 8 16 32 64 128 256];
    log_skip_vals = log2(skip_vals);
    log_skip_vals(1) = 0;
    idx_slice = 5;
    REP = 10000;
elseif L == 8
    skip_vals = [16 32 64 128 256];
    log_skip_vals = log2(skip_vals);
    idx_slice = 9;
    REP = 1000;
end

run_time_avg = zeros(length(skip_vals), length(type_vals));

checkerboard = cell(length(skip_vals), length(type_vals));   % checkerboard - 2
slice = cell(length(skip_vals), length(type_vals));          % slice - 2*L
zerozero = cell(length(skip_vals), length(type_vals));       % zerozero - E, M = 0

mu = cell(1, length(type_vals));
sigma = cell(1, length(type_vals));

for idx_type = 1:length(type_vals)
    type = type_vals(idx_type);
    
    load("Results Average/skip/" + type + "_JDOS_2D_" + L + "L_10E" + log10(REP) + "_nRPS_CPP.mat")
    
    for idx_skip = 1:length(skip_vals)
        skip = skip_vals(idx_skip);
        
        for run = 1:1000
            checkerboard{idx_skip, idx_type}(run) = JDOS{run, idx_skip}(N_SPINS + 1, N_SPINS / 2 + 1);
            slice{idx_skip, idx_type}(run) = JDOS{run, idx_skip}(idx_slice, N_SPINS / 2 + 1);
            zerozero{idx_skip, idx_type}(run) = JDOS{run, idx_skip}(N_SPINS / 2 + 1, N_SPINS / 2 + 1);
        end
        
        run_time_avg(idx_skip, idx_type) = mean(run_time(:, idx_skip));
        
        pd_chkbrd = fitdist(checkerboard{idx_skip, idx_type}', 'normal');
        pd_slice = fitdist(slice{idx_skip, idx_type}', 'normal');
        pd_zerozero = fitdist(zerozero{idx_skip, idx_type}', 'normal');
        
        mu{idx_type}(idx_skip, 1) = pd_chkbrd.mu;
        mu{idx_type}(idx_skip, 2) = pd_slice.mu;
        mu{idx_type}(idx_skip, 3) = pd_zerozero.mu;
        
        sigma{idx_type}(idx_skip, 1) = pd_chkbrd.sigma;
        sigma{idx_type}(idx_skip, 2) = pd_slice.sigma;
        sigma{idx_type}(idx_skip, 3) = pd_zerozero.sigma;
    end
end

if L == 4
    load("JDOS_exact_L4.mat")
    
    mean_error = zeros(length(type_vals), length(skip_vals));
    
    for idx_type = 1:length(type_vals)
        type = type_vals(idx_type);
        
        load("Results Average/skip/" + type + "_JDOS_2D_" + L + "L_10E" + log10(REP) + "_nRPS_CPP.mat")
        
        JDOS_error = cell(1000, length(skip_vals));
        error = zeros(1000, length(skip_vals));
        
        for idx_skip = 1:length(skip_vals)
            skip = skip_vals(idx_skip);
            
            for run = 1:1000
                JDOS_error{run, idx_skip} = abs(JDOS{run, idx_skip}(:, 1:N_SPINS / 2 + 1) - JDOS_exact(:, 1:N_SPINS / 2 + 1)) ./ JDOS_exact(:, 1:N_SPINS / 2 + 1);
                error(run, idx_skip) = sum(sum(JDOS_error{run, idx_skip}(~isnan(JDOS_error{run, idx_skip}))));
            end
            
            mean_error(idx_type, idx_skip) = mean(error(:, idx_skip));
        end
    end
    
    figure(1)
%     subplot(1, 2, 1)
%     plot(skip_vals, mean_error(1, :), '-o')
%     hold on
%     plot(skip_vals, mean_error(2, :), '-o')
%     legend("Error - Par", "Error - Serial", 'location', "best")
%     title("Error Par/Serial vs skip")
%     xlabel("skip")
%     
%     ax2 = axes('Position', [.2 .55 .25 .2]);
%     hold(ax2, 'on')
%     box on
%     plot(ax2, skip_vals(2:5), mean_error(1, 2:5), '-o')
%     hold on
%     plot(ax2, skip_vals(2:5), mean_error(2, 2:5), '-o')
    
%     subplot(1, 2, 2)
    plot(log_skip_vals, mean_error(1, :), '-o')
    hold on
    plot(log_skip_vals, mean_error(2, :), '-o')
    legend("Error - Par", "Error - Serial", 'location', "best")
    title("Error Par/Serial vs skip")
    xlabel("log2(skip)")
end

figure(length(type_vals))
subplot(1, 2, 1)
plot(skip_vals, run_time_avg(:, 1), '-o')
hold on
plot(skip_vals, run_time_avg(:, 2), '-o')
legend("Run Time - Par", "Run Time - Serial", 'location', "NW")
title("Run Time Par/Serial vs skip")
xlabel("skip")

subplot(1, 2, 2)
plot(log_skip_vals, run_time_avg(:, 1), '-o')
hold on
plot(log_skip_vals, run_time_avg(:, 2), '-o')
legend("Run Time - Par", "Run Time - Serial", 'location', "NW")
title("Run Time Par/Serial vs skip")
xlabel("log2(skip)")

figure(length(type_vals) + 1)
subplot(1, 3, 1)
% [plot_checkerboard, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 1), log10(REP_vals), mu{2}(:, 1), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_checkerboard(1).YLim = [0 1.5];
% plot_checkerboard(2).YLim = [0 400];
plot(log_skip_vals, mu{1}(:, 1), '-o')
hold on
plot(log_skip_vals, mu{2}(:, 1), '-o')
legend("mu par", "mu serial", 'location', "best")
title("mu vs skip - checkerboard")
xlabel("log2(skip)")

subplot(1, 3, 2)
% [plot_slice, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 2), log10(REP_vals), mu{2}(:, 2), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_slice(1).YLim = [0 1.5];
% plot_slice(2).YLim = [0 400];
plot(log_skip_vals, mu{1}(:, 2), '-o')
hold on
plot(log_skip_vals, mu{2}(:, 2), '-o')
legend("mu par", "mu serial", 'location', "best")
title("mu vs skip - slice")
xlabel("log2(skip)")

subplot(1, 3, 3)
% [plot_zerozero, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals),mu{1}(:, 3), log10(REP_vals), mu{2}(:, 3), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_zerozero(1).YLim = [0 1.5];
% plot_zerozero(2).YLim = [0 400];
plot(log_skip_vals, mu{1}(:, 3), '-o')
hold on
plot(log_skip_vals, mu{2}(:, 3), '-o')
legend("mu par", "mu serial", 'location', "best")
title("mu vs skip - zerozero")
xlabel("log2(skip)")

figure(length(type_vals) + 2)
subplot(1, 3, 1)
plot(log_skip_vals, sigma{1}(:, 1), '-o')
hold on
plot(log_skip_vals, sigma{2}(:, 1), '-o')
legend("sigma par", "sigma serial", 'location', "best")
title("sigma vs skip - checkerboard")
xlabel("log2(skip)")

subplot(1, 3, 2)
plot(log_skip_vals, sigma{1}(:, 2), '-o')
hold on
plot(log_skip_vals, sigma{2}(:, 2), '-o')
legend("sigma par", "sigma serial", 'location', "best")
title("sigma vs skip - slice")
xlabel("log2(skip)")

subplot(1, 3, 3)
plot(log_skip_vals, sigma{1}(:, 3), '-o')
hold on
plot(log_skip_vals, sigma{2}(:, 3), '-o')
legend("sigma par", "sigma serial", 'location', "best")
title("sigma vs skip - zerozero")
xlabel("log2(skip)")


for idx_type = 1:length(type_vals)
    type = type_vals(idx_type);
    
    figure(idx_type + length(type_vals) + 2)
    subplot(1,3,1)
    h_chkbrd = histfit(checkerboard{skip_vals == skip_hist, idx_type}, [], 'normal');
    %     dataObjs_X = findobj(h_chkbrd,'-property','XData');
    %     x1_chkbrd(:,idx_REP) = dataObjs_X(1).XData;
    %     dataObjs_Y = findobj(h_chkbrd,'-property','YData');
    %     y1_chkbrd(:,idx_REP) = dataObjs_Y(1).YData;
    title("Checkerbord configuration " + type + " - skip = " + skip_hist)
    
    subplot(1,3,2)
    h_slice = histfit(slice{skip_vals == skip_hist, idx_type}, [], 'normal');
    %     dataObjs_X = findobj(h_slice,'-property','XData');
    %     x1_slice(:,idx_REP) = dataObjs_X(1).XData;
    %     dataObjs_Y = findobj(h_slice,'-property','YData');
    %     y1_slice(:,idx_REP) = dataObjs_Y(1).YData;
    title("Slice configuration " + type + " - skip = " + skip_hist)
    
    subplot(1,3,3)
    h_zerozero = histfit(zerozero{skip_vals == skip_hist, idx_type}, [], 'normal');
    %     dataObjs_X = findobj(h_zerozero,'-property','XData');
    %     x1_zerozero(:,idx_REP) = dataObjs_X(1).XData;
    %     dataObjs_Y = findobj(h_zerozero,'-property','YData');
    %     y1_zerozero(:,idx_REP) = dataObjs_Y(1).YData;
    title("E=0,M=0 configuration " + type + " - skip = " + skip_hist)
end
















