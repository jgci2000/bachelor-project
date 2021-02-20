function [S_vector_WL, E_WL_old, WL_log_DOS_dynamic, hist_WL] = ...
    function_WL_RW_step_no_scan(S_vector_WL, E_WL_old, nnxpos, nnxneg, nnypos, nnyneg, E_list, hist_WL, WL_log_DOS_dynamic, f)
%
% RW STEP WITH DYNAMIC WEIGH FACTOR IN LOOP
%
% choose an up and down spin to flip
%
pos = find(S_vector_WL(:,1) == 1);
flipped_pos = randsample(pos,1);
%
neg = find(S_vector_WL(:,1) == -1);
%
if length(neg) > 1
    flipped_neg = randsample(neg,1);
else
    flipped_neg = neg;
end
%
S_vector_WL(flipped_pos,1) = -1;
%
delta_E = - S_vector_WL(flipped_pos,1) .* ( ...
    S_vector_WL(nnxpos(flipped_pos),1) + ...
    S_vector_WL(nnxneg(flipped_pos),1) + ...
    S_vector_WL(nnypos(flipped_pos),1) + ...
    S_vector_WL(nnyneg(flipped_pos),1)) ; % energy of bonds to NN
%
E_WL_new = E_WL_old + 2*delta_E; % build the energy matrix
%
S_vector_WL(flipped_neg,1) = 1;
%
delta_E = - S_vector_WL(flipped_neg,1) .* ( ...
    S_vector_WL(nnxpos(flipped_neg),1) + ...
    S_vector_WL(nnxneg(flipped_neg),1) + ...
    S_vector_WL(nnypos(flipped_neg),1) + ...
    S_vector_WL(nnyneg(flipped_neg),1)) ; % energy of bonds to NN
%
E_WL_new = E_WL_new + 2*delta_E;
%
if rand < min( [ exp(WL_log_DOS_dynamic(E_list == E_WL_old, 1) - WL_log_DOS_dynamic(E_list == E_WL_new, 1)), 1]) % accept
    %
    E_WL_old = E_WL_new;
    WL_log_DOS_dynamic(E_list == E_WL_new, 1) = WL_log_DOS_dynamic(E_list == E_WL_new, 1) + log(f);
    hist_WL(E_list == E_WL_new, 1) = hist_WL(E_list == E_WL_new, 1) + 1;
    %
else % reject
    %
    S_vector_WL(flipped_pos,1) = +1; % unflip spin
    S_vector_WL(flipped_neg,1) = -1; % unflip spin
    WL_log_DOS_dynamic(E_list == E_WL_old, 1) = WL_log_DOS_dynamic(E_list == E_WL_old, 1) + log(f);
    hist_WL(E_list == E_WL_old, 1) = hist_WL(E_list == E_WL_old, 1) + 1;
    %
end
