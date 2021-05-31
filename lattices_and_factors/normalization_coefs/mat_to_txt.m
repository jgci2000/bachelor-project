clear
close

% A simple script that converts the NN table from .mat to .txt
% João Inácio, Oct. 7, 2020

N_vals = [4, 8, 9, 16, 27, 32, 54, 64];%, 108, 128, 256, 512, 1024];
S_vals = [11, 11, 11, 11, 11, 11, 11, 11];% 10, 11, 11, 8, 2];

for idx = 1:length(N_vals)
    N = N_vals(idx);
    S_max = S_vals(idx);
    
    for S = 2:S_max
        if N == 8 && S == 10
            continue
        end
        
        load("coefficients_" + N + "d" + S + ".mat")

        writematrix(norm_factor, "../txt/coefficients_" + N + "d" + S + ".txt")
%         if (N == 32 && S >= 4) || N > 32
%             writematrix(norm_factor_log, "txt/coefficients_log_" + N + "d" + S + ".txt")
%         end
    end
end
