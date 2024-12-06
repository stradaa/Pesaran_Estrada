function [log_mu, log_si] =  getTF_baselines(lfp,w_step, Fs, params, fav_ch, fav_trials)
%GETTF_BASELINES Summary of this function goes here
%   Detailed explanation goes here

arguments
    lfp                     % ch x t x trials
    w_step double = 20      % ms
    Fs     double = 1000    % Hz
    params struct = struct()
    fav_ch = [10,100]; % keep this at two or this will break plot for now
    fav_trials = [1,35];
end

% Default
defaultParams = struct('TrialRejectionMethod', 'rms', ...
                       'TrialRejectionFactor', 2, ...
                       'TF_Type', 'stft', ... % stft, tfspec, etc
                       'TF_Smoothing', 10);   % Hz

% Add missing fields to the input params
fields = fieldnames(defaultParams);
for i = 1:numel(fields)
    if ~isfield(params, fields{i})
        params.(fields{i}) = defaultParams.(fields{i});
    end
end

%% check lfp input
[ch_n, ~, trial_n] = size(lfp);

%% base trial rejection
if trial_n > 1
    % we take average voltage across channels and screen on trial_k x time
    [badTrials,~,~,~,~] = screenBadECogChannels_alex(sq(nanmean(lfp))', ...
                                       1, ... % plot
                                       params.TrialRejectionFactor, ...
                                       params.TrialRejectionMethod, ...
                                       'Baseline Trial #');
    % clean trials
    lfp(:,:,badTrials) = [];
end

% update trial number
trial_n = size(lfp,3);

bad_ch = all(lfp == 0, 2); % get bad ch indexes
if any(ismember(find(bad_ch),fav_ch))
    disp('One of the plotted base channels is not a good channel. Change selection for nicer plot.')
end

%% main loop
all_spec = cell(trial_n,1);

for k = 1:trial_n
    k_lfp = sq(lfp(:,:,k));     % trial
    k_lfp(isnan(k_lfp)) = 0;

    k_spec = cell(ch_n,1);      % trial spec

    iCh_spec_f = cell(ch_n,1); % iCh_f 
    iCh_spec_ti = cell(ch_n,1); % iCh_ti (these get removed with parfor)

    parfor iCh = 1:ch_n
        switch params.TF_Type
            case 'stft' % short time fourier transform
                [iCh_spec, iCh_f, iCh_ti] = computeSpectrogram(k_lfp(iCh,:), ...
                                                               Fs, ...
                                                               params.TF_Smoothing, ...
                                                               w_step);    
                % save outputs
                k_spec{iCh} = abs(iCh_spec).^2;

                iCh_spec_f{iCh} = iCh_f;
                iCh_spec_ti{iCh} = iCh_ti;

            case 'dmtspec' % direct multitaper spectral estimate
                [~, ~] = dmtspec(k_lfp(iCh,:), ...
                                            params.Tapers, ... % [N,W]
                                            Fs, ...        % sampling rate
                                            params.Fk, ... % freq range
                                            params.Pad, ... 
                                            params.Pval, ...
                                            params.Flag);
            case 'tfspec'
                [iCh_spec, iCh_f, iCh_ti] = tfspec(k_lfp(iCh,:), ...
                                           params.Tapers, ... % [N,W]
                                           Fs, ...        % sampling rate
                                           params.Dn, ...
                                           params.Fk, ... % freq range
                                           params.Pad, ... 
                                           params.Pval, ...
                                           params.Flag, ...
                                           params.Contflag);

                if all(iCh_spec==0, 'all') % case when ch = 0
                    k_spec{iCh} = iCh_spec;
                else
                    k_spec{iCh} = log(iCh_spec);
                end
                
                iCh_spec_f{iCh} = iCh_f;
                iCh_spec_ti{iCh} = iCh_ti;
        end
    end

    % format
    spec2 = zeros(length(iCh_spec_ti{1}),length(iCh_spec_f{1}),ch_n); % time x freq x channel
    for i = 1:ch_n
        spec2(:,:,i) = k_spec{i};
    end

    % save trial
    all_spec{k} = spec2;
end

%% Format
% format
spec_f = iCh_spec_f{1};
spec_ti = iCh_spec_ti{1};
combined = zeros(length(spec_ti),length(spec_f),ch_n, trial_n);
for trial = 1:trial_n
    for i = 1:ch_n
        if ~ismember(i, find(bad_ch))
            temp = all_spec{trial}(:,:,i);
            combined(:,:,i,trial) = temp;
        end
    end
end

%% get mean and std across trials
combined_xtime = sq(nanmean(combined,1)); % across time
xtrials_mu = nanmean(combined_xtime,3);
xtrials_sigma = std(combined_xtime,[], 3);

%% output
log_mu = xtrials_mu;
log_si = xtrials_sigma;

%% Check on results
b_check = 1;
if b_check
    figure();t = tiledlayout(5,2);
    title(t, 'Sanity Check: Combined spectrogram for baselines')
    for i = 1:2
        nexttile(i) % plot original
        temp = sq(lfp(fav_ch(i),:,fav_trials(i)));
        temp_x = linspace(params.Plot_bn_base(1), params.Plot_bn_base(2), length(temp));
        plot(temp_x, temp);title(sprintf('LFP | Channel %i | Trial %i', fav_ch(i), fav_trials(i)));
        xlabel('Time (ms)');ylabel('Voltage');

        nexttile(i+2); % plot ch spectrogram
        temp2 = combined(:,:,fav_ch(i),fav_trials(i));
        imagesc(spec_ti, spec_f, temp2');axis xy;
        title(sprintf('Log Spectrogram | Channel %i | Trial %i', fav_ch(i), fav_trials(i)));
        xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
        l_freq = get(gca, 'XLim'); bn = params.Plot_bn_base; 
        n = ceil(diff(bn)/50)+1;
        set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
        set(gca, 'XTickLabel', linspace(bn(1), bn(2), n))
        
        nexttile(i+4); % plot average spectrogram for that ch across tim
        temp3 = nanmean(sq(combined(:,:,fav_ch(i),:)),3); % mean across trials
        imagesc(spec_ti, spec_f, temp3');axis xy;
        title(sprintf('Log Spectrogram | Channel %i | Averaged across all %i trials', fav_ch(i), trial_n));
        xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
        set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
        set(gca, 'XTickLabel', linspace(bn(1), bn(2), n))

        nexttile(i+6); % plot the mu across trials
        temp4 = xtrials_mu(:,fav_ch(i));
        temp5 = xtrials_sigma(:,fav_ch(i));
        yyaxis left; plot(spec_f,temp4);ylabel('mu');
        yyaxis right; plot(spec_f, temp5);ylabel('sigma');
        xlabel('Freq (Hz');
        title(['\mu and \sigma averaged across all trials and time | Ch', num2str(fav_ch(i))])
    end
    
    nexttile(i+7, [1,2]); % plot final across all channels
    temp_6 = sq(nanmean(combined_xtime,3));
    imagesc(1:ch_n, spec_f,temp_6);axis xy;colorbar;
    xlabel('Channel #');ylabel('Freq (Hz)');title('Avg. Spectral Enegy by Ch')
    
    set(gcf, 'Visible', 'on', 'Color', 'w')
end

end

