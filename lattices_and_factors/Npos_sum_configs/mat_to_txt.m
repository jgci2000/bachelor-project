clear
close

% A simple script that converts the NN table from .mat to .txt
% João Inácio, Oct. 7, 2020

Npos_vals = [2, 3, 4, 5, 6, 7, 8, 9, 10];
N_vals = [9, 16];%, 64];

for j = 1:length(Npos_vals)
    N_pos = Npos_vals(j);
    
    for idx = 1:length(N_vals)
        N = N_vals(idx);
        
        load("Npos_sum_configs_Npos" + N_pos + "_N_atm" + N + ".mat")
        
        N_poss = zeros(1, length(Npos_sum_configs));
        for i = 1:length(Npos_sum_configs)
            N_poss(i) = numel(Npos_sum_configs{i}) / N_pos;
        end
        
        write_m = nan(length(Npos_sum_configs), N_pos);

        for i = 1:length(Npos_sum_configs)
            
            if numel(Npos_sum_configs{i}) ~= N_pos
                x = 1;
                for row = 1:numel(Npos_sum_configs{i}) / N_pos
                    for col = 1:N_pos
                        
                        write_m(i, x) = Npos_sum_configs{i}(row, col);
                        x = x + 1;
                    end
                end
            else
                write_m(i, 1:N_pos) = Npos_sum_configs{i};
            end
        end
        
        writematrix(N_poss, "txt/sum_configs_Npos" + N_pos + "_N_atm" + N + ".txt", 'Delimiter', ' ')
        writematrix(write_m, "txt/sum_configs_Npos" + N_pos + "_N_atm" + N + ".txt", 'Delimiter', ' ', 'WriteMode', 'append')
    end
end