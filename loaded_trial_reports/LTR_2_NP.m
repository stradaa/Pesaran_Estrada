function [ch_ignore, Lfp, Lfp_base] = LTR_2_NP(meta, RT, t_factor, show_plot, fs, RM)
%LTR_2 Loaded Trial Report 2
%   This is the NP report given a trial file. Gives an overview with
%   basic cleaning and trial examples from a given center out task.
%   
%   Created for trial types:
%       - delayed_reach
%       - delayed_reach_and_saccade (no saccade info... yet)
%
% Alex Estrada - Nov 14, 2024

arguments
    meta struct
    RT double = 50      % reaction time
    t_factor double = 0.5  % rmoutlier factor 
    show_plot = true   % Creates visible figure (not embedded)
    fs = 2500 % Hz
    RM = 'rms' % channel rejection method
end

%% Impedences
% TODO Not sure how these are done tbh

%% Loading Trials
rec_num = size(meta,2);

% load trials
all_trials = [];
for rec = 1:rec_num
    temp_day = meta(rec).day;
    temp_rec = meta(rec).recs;
    try
        temp_trials = loadTrials(temp_day, temp_rec);
        fprintf('Day: %s, Rec: %s - Trials Loaded\n', temp_day, temp_rec);
    catch ME
        fprintf('Error loading trial data')
        return
    end

    % out
    if rec == 1
        all_trials = temp_trials;
    else
        all_trials = [all_trials, temp_trials]; % ignore bad concatenation
    end
end

% aligned to label
TrialType = {all_trials.PyTaskType};
myInd = ismember(TrialType,meta(1).TrialType);
MyTrials = all_trials(myInd);

% filter by RT
trial_RT = [MyTrials.ReachStart] - [MyTrials.Go]; % not using disGo
MyTrials = MyTrials(trial_RT>RT);

%% Load LFP
% check alignment possbility to label
fld_name = meta(1).align_to;
fld_name_updated = checkFieldForNaNs(all_trials, fld_name);

% load
elecs = 1:384; % hardcoded for now

[Lfp] = trialNPLfp(MyTrials,meta(1).tower, elecs, meta(1).np,fld_name_updated,meta(1).bn);
[Lfp_base] = trialNPLfp(MyTrials,meta(1).tower, elecs, meta(1).np,meta(1).base_fld,meta(1).bn_base);

% reformat
Lfp = permute(Lfp, [2,3,1]); 
Lfp_base = permute(Lfp_base, [2,3,1]);

%% Bad Channel Screening
lfp_mean = sq(mean(Lfp,3));
if isempty(meta(1).ch_ignore)
    [badCh, L, U, C, sd] = screenBadECogChannels_alex(lfp_mean, 1, t_factor, RM, 'Channel');
    ch_ignore = find(badCh);    % channels to ignore by id
else
    [~, L, U, C, sd] = screenBadECogChannels_alex(lfp_mean, 0, t_factor);
    ch_ignore = meta(1).ch_ignore;
end

%% PLOT
recsValues = {meta.recs}; recsString = strjoin(recsValues, ', ');
fig = figure();t = tiledlayout(3,3);
title(t, sprintf('Day: %s, Recs: %s\n Total Reaches: %i. Aligned to: %s', ...
                meta(1).day, ...
                recsString, ...
                length(MyTrials), ...
                meta(1).align_to));

% imp
% nexttile; [numRows,numCols] = size(impGridLayout);
% x = repmat(1:numCols,numCols,1); y = repmat(1:numRows,numRows,1)';
% imagesc(impGridLayout); clim([1 100]); colorbar; hold on;
% temp=num2cell(gridLayout); textString= cellfun(@num2str,temp,'UniformOutput',false);
% text(x(:), y(:), textString, 'HorizontalAlignment', 'Center');hold off;axis('square')
% set(gca, 'XTickLabel', [], 'YTickLabel', []);title('Imp Vals Kohms')

% traces
nexttile(5,[1,2]);
fav_ch = [50,100,200,300];
signal = sq(Lfp(fav_ch,:,:));
avg_signal = sq(nanmean(signal,3));
plot(avg_signal');
legend(split(num2str(fav_ch)));
hold on;
bn_array = linspace(meta(1).bn(1), meta(1).bn(2), length(signal));
[~, aligned_idx] = min(abs(bn_array));
xline(aligned_idx, 'r--',{meta(1).align_to});
hold off;
l_freq = get(gca, 'XLim');
n = ceil(diff(meta(1).bn)/100)+1;
set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
set(gca, 'XTickLabel', linspace(meta(1).bn(1), meta(1).bn(2), n))
title('Mean Ch LFP across reaches');xlabel('Time (ms)');ylabel('Voltage');


% traces bases
nexttile(8,[1,2]);
signal_base = sq(Lfp_base(fav_ch,:,:));
avg_signal = sq(nanmean(signal_base,3));
plot(avg_signal');
legend(split(num2str(fav_ch)));
l_freq = get(gca, 'XLim'); n = ceil(diff(meta(1).bn_base)/100)+1;
set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
set(gca, 'XTickLabel', linspace(meta(1).bn_base(1), meta(1).bn_base(2), n))
title('Mean Ch LFP across baselines')
xlabel('Time (ms)');ylabel('Voltage');

% bad ch
nexttile(4); plot(sd, '.'); hold on
plot(ch_ignore, sd(ch_ignore), 'r*')
yline([L U C], ":", ["Lower Threshold", "Upper Threshold", "Center Value"])
xlabel('Channel'); ylabel(RM)
title('Channel Distribution & First Pass');
hold off;

% mean LFP
temp = sq(mean(Lfp,3));
temp(ch_ignore,:) = nan;

nexttile(1);
histogram(temp);xlabel('Mean activity LFP');ylabel('Count');

nexttile(2,[1,2])
imagesc(temp); xlabel('Time (ms)'); ylabel('Channel');colorbar;
title('Original LFP: Avg X-Trials (badCh removed)')
bn = meta(1).bn; n = ceil(diff(bn)/100)+1;l_freq = get(gca, 'XLim');
x_ticks = linspace(l_freq(1), l_freq(2), n);x_ticks_label = linspace(bn(1), bn(2), n);
set(gca, 'XTick', x_ticks);set(gca, 'XTickLabel', x_ticks_label)
[~, idx_closest_to_zero] = min(abs(x_ticks_label)); % Index closest to 0
x_closest_to_zero = l_freq(1) + (idx_closest_to_zero - 1) * diff(l_freq) / (n - 1);
% Add the vertical line at the closest point to 0
xline(x_closest_to_zero, 'r--','Go','LineWidth',1.5);

if show_plot
    set(fig, 'Color', 'w', 'Visible', 'on')
end


%% Helper functions

    function fld_name = checkFieldForNaNs(all_trials, fld_name)
        % Check if the field exists in the structure
        if isfield(all_trials, fld_name)
            % Get the field values
            fieldValues = [all_trials.(fld_name)];
            nanCounts = sum(isnan(fieldValues));
            
            % Alert if there are any NaN values
            if nanCounts ~= 0
                fprintf('Field "%s" contains %s NaN or missing values\n', ...
                    fld_name, num2str(nanCounts));
                
                % Prompt the user to either enter a new field name or continue
                choice = input('Do you want to enter a new field name? (y/n): ', 's');
                
                if strcmpi(choice, 'y')
                    % Get the new field name from the user
                    fld_name = input('Please enter the new field name: ', 's');
                    % Recursive call to check the new field
                    checkFieldForNaNs(all_trials, fld_name);
                else
                    fprintf('Continuing with field "%s" despite NaN values.\n', fld_name);
                end
            else
                fprintf('Field "%s" does not contain any NaN values.\n', fld_name);
            end
        else
            fprintf('Field "%s" is missing in the structure.\n', fld_name);
        end
    end

end

