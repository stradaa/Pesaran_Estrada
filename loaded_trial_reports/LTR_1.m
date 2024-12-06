function LTR_1(meta_file, show_plot, bin_count)
%LTR_1 Loaded Trial Report 1
%   This is the most basic form of report given a trial file. Gives a
%   general plot that informs analysis given behavior file.
%   
%   Created for trial types:
%       - delayed_reach
%       - delayed_reach_and_saccade (no saccade info... yet)
%
% Alex Estrada - Nov 11, 2024

arguments
    meta_file struct
    show_plot = false   % Creates visible figure (not embedded)
    bin_count = 20      % Default bin count for histogram
end

rec_num = length(meta_file);
fld_names = fieldnames(meta_file);
if ~any(ismember(fld_names, 'day')) || ~any(ismember(fld_names, 'recs'))
    disp('Fieldnames of meta file do not match correct format')
    return
end

%% Load trials
all_trials = [];
for rec = 1:rec_num
    temp_day = meta_file(rec).day;
    temp_rec = meta_file(rec).recs;
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

%% RT
fld_names = fieldnames(all_trials); % reset fieldnames 
requiredFields = {'ReachStart', 'disGo', 'disTargsOn', 'StartAq'};

% check if there's an error
flags = zeros(1,4);
for i = 1:length(requiredFields)
    % Check if missing
    fieldName = requiredFields{i};
    if isfield(all_trials, fieldName)
        % Get the field values
        fieldValues = [all_trials.(fieldName)];
        nanCounts = sum(isnan(fieldValues));
        
        % Alert if there are any NaN values
        if nanCounts ~= 0
            flags(i) = 1;
            fprintf('Field "%s" contains %s NaN or missing values\n', ...
                fieldName, num2str(nanCounts));
        end
    else
        fprintf('Field "%s" is missing in the structure.\n', fieldName);
    end
end

% specify trial type
TrialType = {all_trials.PyTaskType};
myInd = ismember(TrialType,meta_file(1).TrialType);
myTrials = all_trials(myInd);

% get times for plot
if flags(2)
    RT = [myTrials.ReachStart] - [myTrials.Go];
    disp('Using "Go" instead...')
else
    RT = [myTrials.ReachStart] - [myTrials.disGo];
end

reach_duration = [myTrials.ReachStop] - [myTrials.ReachStart];

if flags(3)
    aq_to_targon = [myTrials.TargsOn] - [myTrials.StartAq];
    disp('Using "Targs" instead...')
else
    aq_to_targon = [myTrials.disTargsOn] - [myTrials.StartAq];
end

target_counts = [myTrials.Target];

%% plot
% set-up
recsValues = {meta_file.recs}; recsString = strjoin(recsValues, ', ');
fig = figure();t = tiledlayout(2,2); title(t, sprintf('Day: %s, Recs: %s\n Total Reaches: %i', meta_file(1).day, recsString, length(myTrials)));
set(fig, 'Color', 'w');

% tiled
nexttile;histogram(reach_duration, bin_count);title('Reach Stop - Reach Start');xlabel('Time (ms)');ylabel('Count');
nexttile;histogram(RT, bin_count);title('ReachStart - Go');xlabel('Time (ms)');ylabel('Count');
nexttile;histogram(aq_to_targon, bin_count);title('Target Onset - Start Aq.');xlabel('Time (ms)');ylabel('Count');
nexttile;histogram(target_counts);title('Target Counts');xlabel('Target Label');ylabel('Counts');

% show
if show_plot
    set(fig, 'Visible', 'on')
end


