clear
close all
%
% v1 - works, should be improved in <E> and <C> calculations
% v2 - correcting <E> and <C> calculations | works, should still be improved
%
L = 4;
J = 1;
T(:,1) = 8 : -0.1 : 0.1; %
H(:,1) = 0;%1:-0.05:0; %1:-0.5:0;
%
% REP = 1E4; % when to check for histogram flatness
% p = 0.9; % flatness criteria
% log_f_start = 1E0; % starting log value of f
% log_f_min = 1E-8; % target low value of f
%
% END OF USER INPUT
%
N_atm = L^3;
NN = 6;
%
% wspace_filename = ['JDOS_WL_spinS_Ising_Npos', int2str(Npos), '_2D_SS_L', int2str(L), '_p_', num2str(p), '_log_f_start_1E', int2str(log10(log_f_start)), '_log_f_min_1E', int2str(log10(log_f_min)), '_REP_1E', int2str(log10(REP)), '.mat'];
JDOS_filename = 'JDOS_nRPS_Ising_v1_3D_SC_L4_REP_1E4_skip_1.mat';
%
eval(['load ', JDOS_filename])
%
% NORMALIZATION OF MAGNETIZATION AND ENERGY (imposes unitary length spin vector)
%
M_list = M_list ./ max(M_list) .* N_atm;
E_list = E_list ./ max(E_list) .* 1/2 .* N_atm .* NN; 

%
% FREE ENERGY CALCULATIONS
%
disp('Free Energy calculations...')
%
% FREE ENERGY CALCULATIONS NO MAGVOL
%
Z = zeros(length(H), length(T));
Z_M = zeros(length(H), length(T), length(M_list));
F_temp = nan(length(H), length(T), length(M_list));
%
for q = 1:length(M_list(:,1))
    %
    hits = find(JDOS_nRPS(:,q) ~= 0);
    %
    for k = 1:length(hits)
        %
        Z_M(:,:,q) = Z_M(:,:,q) + JDOS_nRPS(hits(k),q).*exp((E_list(hits(k)) - M_list(q).*H)*(-1./T'));
        %
    end
    %
    Z = Z + Z_M(:,:,q);
    F_temp(:,:,q) = repmat(-T', [length(H),1]).*(log(Z_M(:,:,q))) - M_list(q).*repmat(H,[1,length(T)]); %EDITED
    %
end
%
F = permute(F_temp,[3,2,1]);
%
disp('Free Energy calculations completed')
%
% CALCULATION OF AVERAGE THERMODYNAMIC VALUES
%
avg_M = nan(length(H),length(T));
avg_abs_M = nan(length(H),length(T));
avg_M2 = nan(length(H),length(T));
%
for i=1:length(T)
    %
    for j=1:length(H)
        %
        avg_M(j,i) = 0;
        avg_abs_M(j,i) = 0;
        avg_M2(j,i) = 0;
        %
        for q=1:length(M_list(:,1))
            %
            avg_M(j,i) = avg_M(j,i) + (M_list(q)*Z_M(j,i,q))/Z(j,i);
            avg_abs_M(j,i) = avg_abs_M(j,i) + (abs(M_list(q))*Z_M(j,i,q))/Z(j,i);
            avg_M2(j,i) = avg_M2(j,i) + ((M_list(q)^2)*Z_M(j,i,q))/Z(j,i);
        end
        %
    end
end
%
% FIND MINIMUM OF FREE ENERGY
%
temp = nan(length(M_list), 1);
M_minF = nan(length(H), length(T));
minF = nan(length(H), length(T));
%
for i = 1:length(T)
    for j = 1:length(H)
        for q = 1:length(M_list)
            temp(q,1) = F(q,i,j);
        end
        [u,v] = min(temp(:,1));
        M_minF(j,i) = abs(M_list(v));
        minF(j,i) = u;
        clear v
        %
    end
end
%
M2_minF = M_minF.^2;
%
% DATA NORMALIZATION
%
avg_abs_M = avg_abs_M/(N_atm);
avg_M = avg_M/(N_atm);
avg_M2 = avg_M2/(N_atm^2);
M_minF = M_minF/(N_atm);
M2_minF = M2_minF/(N_atm^2);
minF = minF/(N_atm);
F = F/(N_atm);
%
%
% END OF CALCULATIONS
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% PLOTS FOR T DEPENDENCE
%
if length(H)==1
    %
    figure('Name','Initial plots','NumberTitle','off')
    subplot(2,3,1)
    plot(T, avg_abs_M, '.-'),xlabel('T'), ylabel('<|M|>')
    subplot(2,3,2)
    plot(M_list, F, '.-'),xlabel('M'), ylabel('F per spin')
    subplot(2,3,3)
    plot(T, M_minF, '.-'),xlabel('T'), ylabel('M min F')
    subplot(2,3,4)
    plot(E_list, JDOS_nRPS(:,(length(M_list)+1)/2), '.'), xlabel('E'), ylabel('M=0 DOS')
    subplot(2,3,5)
    plot(E_list, log(JDOS_nRPS(:,(length(M_list)+1)/2)), '.'), xlabel('E'), ylabel('M=0 log DOS')
    %
%     if strcmp(dimension,'2D') && strcmp(lattice,'SS') && strcmp(model,'Ising')
%         %
%         hold on
%         %
%         plot(-1/2*N_atm*NN+4*L,log(2*L),'r.')
%         % slice energy and DOS for Ising 2D SS M=0
%         %
%         plot(2*L^2,log(2),'r.')
%         % checkerboard energy and DOS for Ising 2D SS M=0
%         %
%         hold off
%     end
    %
    
    %
    %
    %             Figure('Name','Average thermodynamic properties (Field-independent)','NumberTitle','oFF')
    %             subplot(1,2,1)
    %             plot(T,avg_abs_M,'k.-',T,real(sqrt(avg_M2-avg_M.^2)),'r.-'),xlabel('T'),ylabel('<|M|>')
    %             subplot(1,2,2)
    %             plot(T,avg_M,'.-'),xlabel('T'),ylabel('<M>')
    %
end
%
% PLOTS FOR H DEPENDENCE
%
if length(T)==1
    %
    figure('Name','Average thermodynamic properties (temperature-independent)','NumberTitle','off')
    subplot(1,2,1)
    plot(H,avg_abs_M,'.-'),xlabel('H'),ylabel('<|M|>')
    subplot(1,2,2)
    plot(H,avg_M,'.-'),xlabel('H'),ylabel('<M>')
    %
end
%
% PLOTS FOR (H,T) DEPENDENCE
%
% UNIT RESCALING
%
%
if length(H)>1 && length(T)>1
    %
    Hsob_avgM=nan(length(H),length(T));
    Hsob_M_minF=nan(length(H),length(T));
    %
    for i=1:length(T)
        %
        Hsob_avgM(1:1:length(H),i)=H(:)./(avg_M(1:1:length(H),i));
        Hsob_M_minF(1:1:length(H),i)=H(:)./(M_minF(1:1:length(H),i));
    end
    %
    % MAGNETOCALORIC EFFECT FROM MAXWELL RELATION
    %
    int=nan(length(H)-1,length(T)-1);
%     T_dS=nan(length(T)-1,1);
    dS_area_avgM=nan(length(H)-1,length(T)-1);
    avgM_dS=nan(length(H)-1,length(T)-1);
    avgM2_dS=nan(length(H)-1,length(T)-1);
    %
    int_minF=nan(length(H)-1,length(T)-1);
    dS_area_M_minF=nan(length(H)-1,length(T)-1);
    M_minF_dS=nan(length(H)-1,length(T)-1);
    M2_minF_dS=nan(length(H)-1,length(T)-1);
    %
    i=1;
    while i <= length(T)
        j=length(H)-1;
        while j >= 1
            int(j,i)=trapz(H(length(H):-1:j),avg_M(length(H):-1:j,i));
            int_minF(j,i)=trapz(H(length(H):-1:j),M_minF(length(H):-1:j,i));
            j=j-1;
        end
        i=i+1;
    end
    %
    k=1; % passo em T(i) para calcular areas
    %
    for i=1:(length(T)-k)
%         T_dS(i,1)=(T(i)+T(i+k))/2;
        dS_area_avgM(:,i)=-(int(:,i+k)-int(:,i))/(T(i+k)-T(i));% * 1E-4;
        avgM_dS(1:1:(length(H)-1),i)=(avg_M(1:1:(length(H)-1),i+k)+avg_M(1:1:(length(H)-1),i))./2;
        avgM2_dS(1:1:(length(H)-1),i)=avgM_dS(1:1:(length(H)-1),i).^2;
        %
        dS_area_M_minF(:,i)=-(int_minF(:,i+k)-int_minF(:,i))/(T(i+k)-T(i));% * 1E-4;
        M_minF_dS(1:1:(length(H)-1),i)=(M_minF(1:1:(length(H)-1),i+k)+M_minF(1:1:(length(H)-1),i))./2;
        M2_minF_dS(1:1:(length(H)-1),i)=M_minF_dS(1:1:(length(H)-1),i).^2;
        %
    end
    %

    %
    figure('Name','Average thermodynamic properties','NumberTitle','off')
    subplot(2,4,1)
    plot(T,avg_abs_M,'k.-'),xlabel('T'),ylabel('<|M|>')
    subplot(2,4,2)
    plot(T,avg_M2-avg_abs_M.^2,'r.-'),xlabel('T'),ylabel('Average susceptibility')
    subplot(2,4,3)
    plot(T,avg_M,'.-'),xlabel('T'),ylabel('<M>')
    subplot(2,4,4)
    plot(H,avg_abs_M,'.-'),xlabel('H'),ylabel('<|M|>')
    subplot(2,4,5)
    plot(H,avg_M,'.-'),xlabel('H'),ylabel('<M>')
    subplot(2,4,6)
    plot(avg_M2,Hsob_avgM,'.-'),xlabel('<M^2>'),ylabel('H/<M>')
    subplot(2,4,7)
    plot(TdS,dS_area_avgM,'.-'),xlabel('T'),ylabel('dS')
    subplot(2,4,8)
    plot(avgM2_dS,dS_area_avgM,'.-'),xlabel('M2'),ylabel('dS')
    %
    figure('Name','Free energy and minimum properties','NumberTitle','off')
    subplot(2,4,1)
    plot(M_list,F(:,:,length(H)),'.-'),ylabel('zero-field Free energy per spin')
    subplot(2,4,2)
    plot(T,minF,'.-'),xlabel('T'),ylabel('Free energy per spin')
    subplot(2,4,3)
    plot(T,M_minF,'.-'),xlabel('T'),ylabel('MminF')
    subplot(2,4,4)
    plot(H,M_minF,'.-'),xlabel('H'),ylabel('MminF')
    subplot(2,4,5)
    plot(M2_minF,Hsob_M_minF,'.-'),xlabel('<MminF^2>'),ylabel('H/MminF')
    subplot(2,4,6)
    plot(TdS,dS_area_M_minF,'.-'),xlabel('T'),ylabel('dS from MminF')
    subplot(2,4,7)
    plot(M2_minF_dS,dS_area_M_minF,'.-'),xlabel('MminF2'),ylabel('dS from MminF')
    %
end
%
disp('Energy and specific heat calculations...')
%
% E_temp = [];
C = nan(length(H),length(T));
%
% for q = 1:length(M_list(:,1))
%     E_temp = vertcat(E_temp,output_parte1_2{q}(:,2));
% end
% %
% E_temp = unique(E_temp);
% E_temp = E_list;
E_temp = E_list(sum(JDOS_nRPS,2)~=0);
%
parte1_final_E=cell(1,length(E_temp));
%
for b = 1:length(E_temp(:,1))
    %
    counter = 0;
    %
    for q = 1:length(M_list(:,1))
        %
        % for z = 1:length(output_parte1_2{q}(:,2))
        for z = 1:length(E_list(:,1))
            %
            % if output_parte1_2{q}(z,2)==E_temp(b,1)
            if E_list(z) == E_temp(b,1)
                %
                counter = counter + 1;
                % parte1_final_E{b}(counter,1) = output_parte1_2{q}(z,1);
                % parte1_final_E{b}(counter,2) = output_parte1_2{q}(z,4);
                %
                parte1_final_E{b}(counter,1) = M_list(q);
                parte1_final_E{b}(counter,2) = JDOS_nRPS(z,q);
                %
            end
        end
    end
end
%
avg_E=nan(length(H),length(T));
avg_E2=nan(length(H),length(T));
%
for j = 1:length(H)
    %
    for i = 1:length(T)
        %
        avg_E(j,i) = 0;
        avg_E2(j,i) = 0;
        %
        for b = 1:length(E_temp(:,1))
            %
            for k = 1:length(parte1_final_E{b}(:,1)) %Q(E)
                %
                avg_E(j,i) = avg_E(j,i) + ((E_temp(b)-parte1_final_E{b}(k,1)*H(j))*parte1_final_E{b}(k,2).*exp(-1/T(i)*(E_temp(b)-parte1_final_E{b}(k,1)*H(j))))/Z(j,i);
                avg_E2(j,i) = avg_E2(j,i) + (((E_temp(b)-parte1_final_E{b}(k,1)*H(j))^2)*parte1_final_E{b}(k,2).*exp(-1/T(i)*(E_temp(b)-parte1_final_E{b}(k,1)*H(j))))/Z(j,i);
                %
            end
            %
        end
        %
        C(j,i) = 1/(T(i)^2).*(avg_E2(j,i) - avg_E(j,i).^2);
        %
    end
    %
end
%
% DATA NORMALIZATION
%
avg_E = avg_E/N_atm;
avg_E2 = avg_E2/N_atm;
C = C/N_atm;
%
disp('Energy and specific heat calculations completed')
%
figure('Name','Energy and Specific Heat','NumberTitle','off')
subplot(1,3,1)
plot(T, avg_E, '.-'), xlabel('T'), ylabel('<Energy> per spin')
subplot(1,3,2)
plot(T, avg_E2, '.-'), xlabel('T'), ylabel('<Energy^2> per spin')
subplot(1,3,3)
plot(T, C, '.-'), xlabel('T'), ylabel('Average specific heat per spin')
%

