% close all
%
figure(1)
%hist_EM(hist_EM == 0) = nan;
L=4;
N_atm=4*4;
NN=4;

E_list=-32:4:32;
M_list=-16:2:16;

load("JDOS_L4_2D_SS_exact")

subplot(1, 2, 2)

b = bar3(JDOS_exact);
zlabel('JDOS')
view(0, 90) % top
%view(133,33)
axis([0 length(M_list)+1 0 length(E_list)+1])
xticks([1 (N_atm/4+1) (N_atm/2+1) (3*N_atm/4+1) N_atm+1])
yticks([1 ((length(E_list)-1)/4 + 1) ((length(E_list)-1)/2+1) (3*(length(E_list)-1)/4+1) length(E_list)])
%
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'XTickLabel', linspace(M_list(min(xt)), M_list(max(xt)), numel(xt))), xlabel('M')
set(gca, 'YTickLabel', linspace(E_list(min(yt)), E_list(max(yt)), numel(yt))), ylabel('E')
%
h = get(gca,'children');
%
for w = 1:length(h)
    %
    hc = get(h(w),'cdata');
    hz = get(h(w),'zdata');
    %
    for u = 1:(length(hc(:,1))/6)
        %
        if sum(nansum( hc(1+(u-1)*6 : u*6, 1:4) )) == 0
            %
            hc_new = hc;
            hc_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'cdata', hc_new);
            hc = hc_new;
            %
            hz_new = hz;
            hz_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'zdata', hz_new);
            hz = hz_new;
            %
        end
        %
    end
    %
end
%
cbh = colorbar;
cbh.Label.String = 'JDOS';
