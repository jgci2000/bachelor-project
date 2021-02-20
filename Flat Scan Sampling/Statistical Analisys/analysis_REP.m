clear
close

% Error analysis for the CPP results
% João Inácio, 12 Nov., 2020

% Change REP_hist to change which REP value is ploted in the histogram
% plots

L = 8;
N_SPINS = L * L;
skip = N_SPINS;
type_vals = ["Par", "Serial"]; % Serial/Par

REP_hist = 10000;

if L == 4
    REP_vals = [100, 1000, 10000, 100000, 1000000];
    idx_slice = 5;
elseif L == 8
    REP_vals = [1000, 10000];
    idx_slice = 9;
end

run_time_avg = zeros(length(REP_vals), length(type_vals));

checkerboard = cell(length(REP_vals), length(type_vals));   % checkerboard - 2
slice = cell(length(REP_vals), length(type_vals));          % slice - 2*L
zerozero = cell(length(REP_vals), length(type_vals));       % zerozero - E, M = 0

mu = cell(1, length(type_vals));
sigma = cell(1, length(type_vals));

for idx_type = 1:length(type_vals)
    type = type_vals(idx_type);
    
    load("Results Average/REP/" + type + "_JDOS_2D_" + L + "L_skip" + skip + "_nRPS_CPP.mat")
    
    for idx_REP = 1:length(REP_vals)
        REP = REP_vals(idx_REP);
        
        for run = 1:1000
            checkerboard{idx_REP, idx_type}(run) = JDOS{run, idx_REP}(N_SPINS + 1, N_SPINS / 2 + 1);
            slice{idx_REP, idx_type}(run) = JDOS{run, idx_REP}(idx_slice, N_SPINS / 2 + 1);
            zerozero{idx_REP, idx_type}(run) = JDOS{run, idx_REP}(N_SPINS / 2 + 1, N_SPINS / 2 + 1);
        end
        
        run_time_avg(idx_REP, idx_type) = mean(run_time(:, idx_REP));
        
        pd_chkbrd = fitdist(checkerboard{idx_REP, idx_type}', 'normal');
        pd_slice = fitdist(slice{idx_REP, idx_type}', 'normal');
        pd_zerozero = fitdist(zerozero{idx_REP, idx_type}', 'normal');
        
        mu{idx_type}(idx_REP, 1) = pd_chkbrd.mu;
        mu{idx_type}(idx_REP, 2) = pd_slice.mu;
        mu{idx_type}(idx_REP, 3) = pd_zerozero.mu;
        
        sigma{idx_type}(idx_REP, 1) = pd_chkbrd.sigma;
        sigma{idx_type}(idx_REP, 2) = pd_slice.sigma;
        sigma{idx_type}(idx_REP, 3) = pd_zerozero.sigma;
    end
end

if L == 4
    load("JDOS_exact_L4.mat")
    
    mean_error = zeros(length(type_vals), length(REP_vals));
    
    for idx_type = 1:length(type_vals)
        type = type_vals(idx_type);
        
        load("Results Average/REP/" + type + "_JDOS_2D_" + L + "L_skip" + skip + "_nRPS_CPP.mat")
        
        JDOS_error = cell(1000, length(REP_vals));
        error = zeros(1000, length(REP_vals));
        
        for idx_REP = 1:length(REP_vals)
            REP = REP_vals(idx_REP);
            
            for run = 1:1000
                JDOS_error{run, idx_REP} = abs(JDOS{run, idx_REP}(:, 1:N_SPINS / 2 + 1) - JDOS_exact(:, 1:N_SPINS / 2 + 1)) ./ JDOS_exact(:, 1:N_SPINS / 2 + 1);
                error(run, idx_REP) = sum(sum(JDOS_error{run, idx_REP}(~isnan(JDOS_error{run, idx_REP}))));
            end
            
            mean_error(idx_type, idx_REP) = mean(error(:, idx_REP));
        end
        
    end
    
    figure(1)
%     subplot(1, 2, 1)
%     plot(REP_vals, mean_error(1, :), '-o')
%     hold on
%     plot(REP_vals, mean_error(2, :), '-o')
%     legend("Error - Par", "Error - Serial", 'location', "best")
%     title("Error Par/Serial vs REP")
%     xlabel("REP")
%     
%     subplot(1, 2, 2)
    plot(log10(REP_vals), mean_error(1, :), '-o')
    hold on
    plot(log10(REP_vals), mean_error(2, :), '-o')
    legend("Error - Par", "Error - Serial", 'location', "best")
    title("Error Par/Serial vs REP")
    xlabel("log(REP)")
end

figure(length(type_vals))
subplot(1, 2, 1)
plot(REP_vals, run_time_avg(:, 1), '-o')
hold on
plot(REP_vals, run_time_avg(:, 2), '-o')
legend("Run Time - Par", "Run Time - Serial", 'location', "NW")
title("Run Time Par/Serial vs REP")
xlabel("REP")

subplot(1, 2, 2)
plot(log10(REP_vals), run_time_avg(:, 1), '-o')
hold on
plot(log10(REP_vals), run_time_avg(:, 2), '-o')
legend("Run Time - Par", "Run Time - Serial", 'location', "NW")
title("Run Time Par/Serial vs REP")
xlabel("log(REP)")

figure(length(type_vals) + 1)
subplot(1, 3, 1)
% [plot_checkerboard, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 1), log10(REP_vals), mu{2}(:, 1), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_checkerboard(1).YLim = [0 1.5];
% plot_checkerboard(2).YLim = [0 400];
plot(log10(REP_vals), mu{1}(:, 1), '-o')
hold on
plot(log10(REP_vals), mu{2}(:, 1), '-o')
legend("mu par", "mu serial", 'location', "best")
title("mu vs REP - checkerboard")
xlabel("log(REP)")

subplot(1, 3, 2)
% [plot_slice, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals), mu{1}(:, 2), log10(REP_vals), mu{2}(:, 2), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_slice(1).YLim = [0 1.5];
% plot_slice(2).YLim = [0 400];
plot(log10(REP_vals), mu{1}(:, 2), '-o')
hold on
plot(log10(REP_vals), mu{2}(:, 2), '-o')
legend("mu par", "mu serial", 'location', "best")
title("mu vs REP - slice")
xlabel("log(REP)")

subplot(1, 3, 3)
% [plot_zerozero, mu_par_line, mu_serial_line] = plotyy(log10(REP_vals),mu{1}(:, 3), log10(REP_vals), mu{2}(:, 3), 'plot', 'plot');
% mu_par_line.Marker = 'o';
% mu_serial_line.Marker = 'o';
% plot_zerozero(1).YLim = [0 1.5];
% plot_zerozero(2).YLim = [0 400];
plot(log10(REP_vals), mu{1}(:, 3), '-o')
hold on
plot(log10(REP_vals), mu{2}(:, 3), '-o')
legend("mu par", "mu serial", 'location', "best")
title("mu vs REP - zerozero")
xlabel("log(REP)")

figure(length(type_vals) + 2)
subplot(1, 3, 1)
plot(log10(REP_vals), sigma{1}(:, 1), '-o')
hold on
plot(log10(REP_vals), sigma{2}(:, 1), '-o')
legend("sigma par", "sigma serial", 'location', "best")
title("sigma vs REP - checkerboard")
xlabel("log(REP)")

subplot(1, 3, 2)
plot(log10(REP_vals), sigma{1}(:, 2), '-o')
hold on
plot(log10(REP_vals), sigma{2}(:, 2), '-o')
legend("sigma par", "sigma serial", 'location', "best")
title("sigma vs REP - slice")
xlabel("log(REP)")

subplot(1, 3, 3)
plot(log10(REP_vals), sigma{1}(:, 3), '-o')
hold on
plot(log10(REP_vals), sigma{2}(:, 3), '-o')
legend("sigma par", "sigma serial", 'location', "best")
title("sigma vs REP - zerozero")
xlabel("log(REP)")


for idx_type = 1:length(type_vals)
    type = type_vals(idx_type);
    
    figure(idx_type + length(type_vals) + 2)
    subplot(1,3,1)
    h_chkbrd = histfit(checkerboard{REP_vals == REP_hist, idx_type}, [], 'normal');
    %     dataObjs_X = findobj(h_chkbrd,'-property','XData');
    %     x1_chkbrd(:,idx_REP) = dataObjs_X(1).XData;
    %     dataObjs_Y = findobj(h_chkbrd,'-property','YData');
    %     y1_chkbrd(:,idx_REP) = dataObjs_Y(1).YData;
    title("Checkerbord configuration " + type + " - REP = " + REP_hist)
    
    subplot(1,3,2)
    h_slice = histfit(slice{REP_vals == REP_hist, idx_type}, [], 'normal');
    %     dataObjs_X = findobj(h_slice,'-property','XData');
    %     x1_slice(:,idx_REP) = dataObjs_X(1).XData;
    %     dataObjs_Y = findobj(h_slice,'-property','YData');
    %     y1_slice(:,idx_REP) = dataObjs_Y(1).YData;
    title("Slice configuration " + type + " - REP = " + REP_hist)
    
    subplot(1,3,3)
    h_zerozero = histfit(zerozero{REP_vals == REP_hist, idx_type}, [], 'normal');
    %     dataObjs_X = findobj(h_zerozero,'-property','XData');
    %     x1_zerozero(:,idx_REP) = dataObjs_X(1).XData;
    %     dataObjs_Y = findobj(h_zerozero,'-property','YData');
    %     y1_zerozero(:,idx_REP) = dataObjs_Y(1).YData;
    title("E=0,M=0 configuration " + type + " - REP = " + REP_hist)
end
















