function [combined_spec, combined_spec_norm] = getTF_Spec(lfp, lfp_bases, w_step, Fs, params, show_plot)
%GETTF_SPEC Summary of this function goes here
%   Temporary before creating a more generalizable one

arguments
    lfp                     % channel x time x trial
    lfp_bases               % channel x time x trial
    w_step double = 20      % ms
    Fs     double = 1000    % Hz (default ECoG 1000Hz or NP 2500 for cLfp)
    params struct = struct()
    show_plot = 1
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

% check lfp input
[ch_n, ~, trial_n] = size(lfp);

% get baseline
[mu,sigma] = getTF_baselines(lfp_bases,w_step, Fs, params);

%% Main loop
all_spec = cell(trial_n,1);
all_spec_norm = cell(trial_n,1);

for k = 1:trial_n
    k_lfp = sq(lfp(:,:,k));     % trial
    k_lfp(isnan(k_lfp)) = 0;

    k_spec = cell(ch_n,1);      % trial spec
    iCh_spec_f = cell(ch_n,1);  % iCh_f 
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
                [iCh_spec, iCh_f] = dmtspec(k_lfp(iCh,:), ...
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
            
    % normalize by trial bc dealing with 4D plots is a mess
    baseline_mu_expanded = permute(repmat(mu, 1, 1, length(iCh_spec_ti{1})), [3, 1, 2]);
    baseline_sigma_expanded = permute(repmat(sigma, 1, 1, length(iCh_spec_ti{1})), [3, 1, 2]);

    % spec_subMean = spec2 - baseline_mu_expanded; % no division by variance
    zspec_norm = (spec2 - baseline_mu_expanded)./baseline_sigma_expanded;

    % save trial
    all_spec{k} = spec2;
    all_spec_norm{k} = zspec_norm;
end

%% Format
spec_f = iCh_spec_f{1};
spec_ti = iCh_spec_ti{1};
combined_spec = zeros(length(spec_ti),length(spec_f),ch_n, trial_n);
combined_spec_norm = zeros(length(spec_ti),length(spec_f),ch_n, trial_n);

for trial = 1:trial_n
    for i = 1:ch_n
        if ~ismember(i, params.bad_ch)
            temp = all_spec{trial}(:,:,i);
            temp_norm = all_spec_norm{trial}(:,:,i);

            combined_spec(:,:,i,trial) = temp;
            combined_spec_norm(:,:,i,trial) = temp_norm;
        end
    end
end

%% Plot
if show_plot
    %% Figure 1
    fig = figure();t = tiledlayout(3,1);

    tr = 1;
    ch = 50;

    title(t, sprintf('Sanity Check | Reach %i | Ch %i',tr,ch));
    
    time_trace = sq(lfp(ch,:,tr));
    spec_ori = sq(combined_spec(:,:,ch,tr));
    spec_norm = sq(combined_spec_norm(:,:,ch,tr));
    
    nexttile(1); % plot time trace
    bn = params.Plot_bn;
    x_trace = linspace(bn(1), bn(2), length(time_trace));
    [~, zero_idx] = min(abs(x_trace));zero_x = x_trace(zero_idx); % Corresponding x value
    plot(x_trace, time_trace);hold on;xline(zero_x, 'r--','Go','LineWidth',1.5);
    xlabel('Time (ms)');ylabel('Voltage');title('LFP: Time trace')

    nexttile(2);
    imagesc(spec_ti, spec_f, spec_ori');axis xy;title('Original Log Spectrogram');
    xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
    l_freq = get(gca, 'XLim');n = ceil(diff(bn)/100)+1;
    x_ticks = linspace(l_freq(1), l_freq(2), n);x_ticks_label = linspace(bn(1), bn(2), n);
    set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
    [~, idx_closest_to_zero] = min(abs(x_ticks_label)); % Index closest to 0
    x_closest_to_zero = l_freq(1) + (idx_closest_to_zero - 1) * diff(l_freq) / (n - 1);
    % Add the vertical line at the closest point to 0
    xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);

    nexttile(3);
    imagesc(spec_ti, spec_f, spec_norm');axis xy;title('Z-Norm Log Spectrogram');
    xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
    set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
    set(gca, 'XTickLabel', linspace(bn(1), bn(2), n))
    xline(x_closest_to_zero, 'r--','Go', 'LineWidth',1.5);

    % end of figure 1
    set(fig, 'Color', 'w', 'Visible', 'on')
    clear spec_norm spec_ori

    %% Figure 2 (Means etc)
    fig = figure();t = tiledlayout(3,1);title(t, 'Averages');
    
    xtrials_lfp = sq(nanmean(lfp,3));

    xtrials_spec = nanmean(combined_spec,4);
    xtrials_xch_spec = nanmean(xtrials_spec,3);

    xtrials_spec_norm = nanmean(combined_spec_norm,4);
    xtrials_xch_spec_norm = nanmean(xtrials_spec_norm,3);
    
    nexttile(1);
    imagesc(xtrials_lfp);title('Avg. LFP X-trials');
    xlabel('Time (ms)');ylabel('Channel #');
    l_freq = get(gca, 'XLim');n = ceil(diff(bn)/100)+1;
    x_ticks = linspace(l_freq(1), l_freq(2), n);x_ticks_label = linspace(bn(1), bn(2), n);
    set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
    [~, idx_closest_to_zero] = min(abs(x_ticks_label)); % Index closest to 0
    x_closest_to_zero = l_freq(1) + (idx_closest_to_zero - 1) * diff(l_freq) / (n - 1);
    % Add the vertical line at the closest point to 0
    xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);
    
    nexttile(2);
    imagesc(spec_ti, spec_f, xtrials_xch_spec');axis xy;title('Mean X-trials & X-Ch Original Log Spectrogram');
    xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
    l_freq = get(gca, 'XLim');n = ceil(diff(bn)/100)+1;
    x_ticks = linspace(l_freq(1), l_freq(2), n);x_ticks_label = linspace(bn(1), bn(2), n);
    set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
    [~, idx_closest_to_zero] = min(abs(x_ticks_label)); % Index closest to 0
    x_closest_to_zero = l_freq(1) + (idx_closest_to_zero - 1) * diff(l_freq) / (n - 1);
    % Add the vertical line at the closest point to 0
    xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);

    nexttile(3);
    imagesc(spec_ti, spec_f, xtrials_xch_spec_norm');axis xy;title('Mean X-trials & X-Ch Norm Log Spectrogram');
    xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
    set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
    % Add the vertical line at the closest point to 0
    xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);
    
    % ARTIFACT SETTING YLIM
    new_y_range = 1;
    if new_y_range
        new_ylim = [70,150];
        y_idx = (spec_f >= new_ylim(1) & spec_f <= new_ylim(2));
    
        nexttile(2);
        imagesc(spec_ti, spec_f(y_idx), xtrials_xch_spec(:,y_idx)');axis xy;title('Mean X-trials & X-Ch Original Log Spectrogram');
        xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
        set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
        % Add the vertical line at the closest point to 0
        xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);
    
        nexttile(3);
        imagesc(spec_ti, spec_f(y_idx), xtrials_xch_spec_norm(:,y_idx)');axis xy;title('Mean X-trials & X-Ch Norm Log Spectrogram');
        xlabel('Time (ms)');ylabel('Freq (Hz)');colorbar;
        set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
        % Add the vertical line at the closest point to 0
        xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);
    end
    
    set(fig, 'Color', 'w', 'Visible', 'on')
    clear xtrials_lfp xtrials_spec xtrials_spec_norm xtrials_xch_spec xtrials_xch_spec_norm
    
    %% Tuning Curves
    load(['myTrials_',params.day, '.mat'])
    targets = [myTrials.ChosenTarget];

    figure(); title('Tuning Curves');

    % per direction
    for direction = 1:7
        idx = (targets==direction+1);
        
        spec_targ = nanmean(combined_spec_norm(:,:,:,idx),4); % mean by trial
        spec_targ2 = nanmean(spec_targ,3); % mean by channel
        
        if new_y_range
            spec_targ3 = nanmean(spec_targ2(:,y_idx),2); % average by select Freq
        else
            spec_targ3 = nanmean(spec_targ2,2); % average by all freq
        end

        x_trace = linspace(bn(1), bn(2), length(spec_targ3));
        plot(x_trace, spec_targ3);
        hold on;
    end
    legend(split(num2str(2:8)))
    xline(0, 'r--','Go','LineWidth',1.5);
    hold off;xlabel('Time (ms)');ylabel('Norm Power');
    set(gcf, 'Color', 'w')
end

end

