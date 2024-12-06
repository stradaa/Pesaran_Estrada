function [ch_ignore, lfp_reach_ref, lfp_base_ref] = LTR_2_ECoG(meta, RT, t_factor, reref_type, show_plot, fs)
%LTR_2 Loaded Trial Report 2
%   This is the uECoG report given a trial file. Gives an overview with
%   basic cleaning and trial examples from a given center out task.
%   
%   Created for trial types:
%       - delayed_reach
%       - delayed_reach_and_saccade (no saccade info... yet)
%
% Alex Estrada - Nov 14, 2024

arguments
    meta struct
    RT double = 50     % reaction time
    t_factor double = 0.5  % rmoutlier factor 
    reref_type {string, char} = 'CARLA'   % type of re-ref (add more)
    show_plot = true   % Creates visible figure (not embedded)
    fs = 1000 % Hz
end

%% Impedences
[~, ~, gridLayout, impGridLayout] = getImpVals_ECoG_alex(meta);

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

% load reach
Lfp_Reach = trialLfp(MyTrials, ...
                     meta(1).tower, ...
                     1:244, ...
                     1, ...
                     fld_name_updated, ...
                     meta(1).bn, ...
                     meta(1).MONKEYDIR, ...
                     meta(1).lfpType);
% load bases
got_bases = 0;
if isfield(meta(1), 'base_fld')
    Lfp_Base = trialLfp(MyTrials, ...
                        meta(1).tower, ...
                        1:244, ...
                        1, ...
                        meta(1).base_fld, ...
                        meta(1).bn_base, ...
                        meta(1).MONKEYDIR, ...
                        meta(1).lfpType);
    got_bases= 1;
end

%% Bad Channel Screening
lfp_mean = sq(mean(Lfp_Reach,1));
if isempty(meta(1).ch_ignore)
    [badCh, L, U, C, sd] = screenBadECogChannels_alex(lfp_mean, 0, t_factor);
    ch_ignore = find(badCh);    % channels to ignore by id
else
    [~, L, U, C, sd] = screenBadECogChannels_alex(lfp_mean, 0, t_factor);
    ch_ignore = meta(1).ch_ignore;
end

%% Re-referencing
len_reach = size(lfp_mean,2);
time_start = tic; fprintf('Starting re-referncing procedure: %s\n', reref_type);
switch reref_type
    case 'bypcbid'
        lfp_reach_ref = zeros(size(Lfp_Reach));
        for trial = 1:size(Lfp_Reach, 1)
            data_temp = sq(Lfp_Reach(trial, :, :));
            % re_reference by pcb_id
            [pcb_reref, ~] = avgref_bypcbid(data_temp, meta(1).experiment, ch_ignore, 1, 0);
            % save
            lfp_reach_ref(trial, :,:) = pcb_reref;
        end

        if got_bases
            lfp_base_ref = zeros(size(Lfp_Base));
            for trial = 1:size(Lfp_Reach,1)
                data_temp = sq(Lfp_Base(trial,:,:));
                % re_reference by pcb_id
                [pcb_reref, ~] = avgref_bypcbid(data_temp, meta(1).experiment, ch_ignore, 1, 0);
                % save
                lfp_base_ref(trial, :,:) = pcb_reref;
            end
        end

    case 'CARLA'
        tt = linspace(0,len_reach/fs, len_reach);
        data_temp = permute(Lfp_Reach, [2,3,1]);    % reshape to specified format
        lfp_reach_ref = zeros(size(data_temp)); % initialize with format
        % create mask
        ch_mask = true(size(data_temp,1), 1);
        ch_mask(ch_ignore) = false;
        % remove bad channels
        data_temp = data_temp(ch_mask,:,:);
        % CARLA
        [out_temp, CAR, stats] = CARLA(tt, data_temp,fs);
        % repopulate
        lfp_reach_ref(ch_mask,:,:) = out_temp;
        % reshape
        lfp_reach_ref = permute(lfp_reach_ref, [3,1,2]);
end
fprintf("Complete. Time Elapsed: %.2f seconds", toc(time_start))

%% PLOT
recsValues = {meta.recs}; recsString = strjoin(recsValues, ', ');
fig = figure();t = tiledlayout(3,3);
title(t, sprintf('Day: %s, Recs: %s\n Total Reaches: %i. Aligned to: %s', ...
                meta(1).day, ...
                recsString, ...
                length(MyTrials), ...
                meta(1).align_to));

% imp
nexttile; [numRows,numCols] = size(impGridLayout);
x = repmat(1:numCols,numCols,1); y = repmat(1:numRows,numRows,1)';
imagesc(impGridLayout); clim([1 100]); colorbar; hold on;
temp=num2cell(gridLayout); textString= cellfun(@num2str,temp,'UniformOutput',false);
text(x(:), y(:), textString, 'HorizontalAlignment', 'Center');hold off;axis('square')
set(gca, 'XTickLabel', [], 'YTickLabel', []);title('Imp Vals Kohms')

% bad ch
nexttile(4); plot(sd, '.'); hold on
plot(ch_ignore, sd(ch_ignore), 'r*')
yline([L U C], ":", ["Lower Threshold", "Upper Threshold", "Center Value"])
xlabel('Channel'); ylabel('')
title('RMS: Distribution per Channel');

% avg original vs reref
nexttile(2, [1,2]);
imagesc(sq(mean(Lfp_Reach,1))); xlabel('Time (ms)'); ylabel('Channel'); 
title('Original LFP (averaged across trials)')
l_freq = get(gca, 'XLim');
n = ceil(diff(meta(1).bn)/100)+1;
set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
set(gca, 'XTickLabel', linspace(meta(1).bn(1), meta(1).bn(2), n))
colorbar;

nexttile(5, [1,2]);
imagesc(sq(mean(lfp_reach_ref,1))); xlabel('Time (ms)'); ylabel('Channel'); 
title('Re-Referenced LFP by PCB id (averaged across trials)');
l_freq = get(gca, 'XLim');
n = ceil(diff(meta(1).bn)/100)+1;
set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
set(gca, 'XTickLabel', linspace(meta(1).bn(1), meta(1).bn(2), n))
colorbar;

if got_bases
    nexttile(8, [1,2])
    fav_ch = [35, 150, 200];
    signal_base = sq(lfp_base_ref(:,fav_ch,:));
    avg_signal = sq(nanmean(signal_base));
    plot(avg_signal');
    legend(split(num2str(fav_ch)));
    l_freq = get(gca, 'XLim');
    n = ceil(diff(meta(1).bn_base)/100)+1;
    set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
    set(gca, 'XTickLabel', linspace(meta(1).bn_base(1), meta(1).bn_base(2), n))
    title('Mean Ch LFP across baselines')
    xlabel('Time (ms)');ylabel('Voltage');
end

if show_plot
    set(fig, 'Color', 'w')
    set(fig, 'Visible', 'on')
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

