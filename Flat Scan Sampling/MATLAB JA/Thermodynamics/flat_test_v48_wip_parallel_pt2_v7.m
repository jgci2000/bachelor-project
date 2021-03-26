clear
close all
%
% v1 - first version
% v2 - now correct, from high T to low T until sign change in dif
% v3 - consider all M > 0 values
% v4 - also distribution of chkbrd and slice
% v5 - also distribution of 0,0
% v6 - also Tc of average JDOS
% v7 - optimization of Tc calc, Tc of exact JDOS for L = 4
%
load workspace_ft_v48_parallel_L8R1E4x1E3.mat
%
L = 8;
q_M0 = L^2/2 + 1;
%
T(:,1) = 2.5 : 0.001 : 3.3;
%
chkbrd_all = zeros(floor(REP/slice),1);
slice_all = zeros(floor(REP/slice),1);
zerozero_all = zeros(floor(REP/slice),1);
Tc_all = zeros(floor(REP/slice),1);
%
tic
%
Z_Mplus_temp = nan(q_M0-1, 1);
%
for BIG_K = 1:floor(REP/slice)
    %
    chkbrd_all(BIG_K,1) = JDOS_aprox{BIG_K}( E_list == 1/2 * N_atm * NN, N_atm/2 + 1);
    slice_all(BIG_K,1) = JDOS_aprox{BIG_K}( E_list == -1/2 * NN * L^2 + 2 * 2 * L, N_atm/2 + 1);
    zerozero_all(BIG_K,1) = JDOS_aprox{BIG_K}( E_list == 0, N_atm/2 + 1);
    %
    for  j = length(T):-1:1
        %
        Z_M0_temp = sum(JDOS_aprox{BIG_K}(:,q_M0) .* exp(-E_list/T(j)));
        %
        for q = 1:(q_M0-1)
            %
            Z_Mplus_temp(q,1) = sum(JDOS_aprox{BIG_K}(:,q) .* exp(-E_list/T(j)));
            %
        end
        %
        dif_temp = Z_Mplus_temp - Z_M0_temp;
        %
        if any(dif_temp(:,1) > 0 )
            %
            break
            %
        end
        %
    end
    %
    Tc_all(BIG_K,1) = (T(j)+T(j-1))/2;
    %
end
%
toc
%
figure(1)
subplot(1,3,1)
h_chkbrd = histfit(chkbrd_all, [], 'normal');
dataObjs_X = findobj(h_chkbrd,'-property','XData');
x1_chkbrd(:,1) = dataObjs_X(1).XData;
dataObjs_Y = findobj(h_chkbrd,'-property','YData');
y1_chkbrd(:,1) = dataObjs_Y(1).YData;
%
subplot(1,3,2)
h_slice = histfit(slice_all, [], 'normal');
dataObjs_X = findobj(h_slice,'-property','XData');
x1_slice(:,1) = dataObjs_X(1).XData;
dataObjs_Y = findobj(h_slice,'-property','YData');
y1_slice(:,1) = dataObjs_Y(1).YData;
%
subplot(1,3,3)
h_zerozero = histfit(zerozero_all, [], 'normal');
dataObjs_X = findobj(h_zerozero,'-property','XData');
x1_zerozero(:,1) = dataObjs_X(1).XData;
dataObjs_Y = findobj(h_zerozero,'-property','YData');
y1_zerozero(:,1) = dataObjs_Y(1).YData;
%
pd_chkbrd = fitdist(chkbrd_all, 'normal');
pd_slice = fitdist(slice_all, 'normal');
pd_zerozero = fitdist(zerozero_all, 'normal');
%
figure(2)
h = histfit(Tc_all, [], 'normal');
pd_Tc = fitdist(Tc_all, 'normal');
%
dataObjs_X = findobj(h,'-property','XData');
x1_Tc(:,1) = dataObjs_X(1).XData;
%
dataObjs_Y = findobj(h,'-property','YData');
y1_Tc(:,1) = dataObjs_Y(1).YData;
%
% Tc of average JDOS
%
for  j = length(T):-1:1
    %
    Z_M0_temp_avg = sum(avg_JDOS_aprox(:,q_M0) .* exp(-E_list/T(j)));
    %
    for q = 1:(q_M0-1)
        %
        Z_Mplus_temp_avg(q,1) = sum(avg_JDOS_aprox(:,q) .* exp(-E_list/T(j)));
        %
    end
    %
    dif_temp_avg = Z_Mplus_temp_avg - Z_M0_temp_avg;
    %
    if any(dif_temp_avg(:,1) > 0 )
        %
        break
        %
    end
    %
end
%
Tc_avg = (T(j)+T(j-1))/2;
%
% Tc of exact JDOS for L = 4
%
if L == 4
    %
    load JDOS_L4_2D_SS_exact.mat
    %
    for  j = length(T):-1:1
        %
        Z_M0_temp_exact = sum(JDOS_exact(:,q_M0) .* exp(-E_list/T(j)));
        %
        for q = 1:(q_M0-1)
            %
            Z_Mplus_temp_exact(q,1) = sum(JDOS_exact(:,q) .* exp(-E_list/T(j)));
            %
        end
        %
        dif_temp_exact = Z_Mplus_temp_exact - Z_M0_temp_exact;
        %
        if any(dif_temp_exact(:,1) > 0 )
            %
            break
            %
        end
        %
    end
    %
    Tc_exact = (T(j)+T(j-1))/2;
    %
end
