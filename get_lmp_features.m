function [outputArg1,outputArg2] = get_lmp_features(inputArg1,inputArg2)
%GET_LMP_FEATURES Summary of this function goes here
%   Detailed explanation goes here
arguments
    Data_Reach
    Data_Base
    baseline
    params
    to_plot
end


% baseline
base_permuted = permute(Data_Base, [2, 3, 1]); % channels x time x trials
base_reshaped = reshape(base_permuted, 244, []); % channels x (time * trials)

% lmp extraction
lmp_base = filtfilt(n, Wn, base_reshaped')';
% normalize values 
[z,mu,sigma] = zscore(lmp_base');

% filter
for i = 1:size(Data_Reach,1)
    reach = sq(Data_Reach(i,:,:));

    % applying filter
    lmp = filtfilt(n, Wn, reach')';
    
    % normalize
    [z,mu,sigma] = zscore(lmp_base');
    lmp_norm = (lmp - mu') ./ sigma';

    % save
    lmp_filtfilt(i,:,:) = lmp_norm;

    if to_plot
        plot_normalization(reach, lmp, lmp_norm, i, t_used);
    end
end
end

