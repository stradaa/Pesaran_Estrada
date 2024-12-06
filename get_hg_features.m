function all_formatted = get_hg_features(data, spec_params, base, base_type, to_plot, x_tics_n)
%GET_HG_FEATURES Summary of this function goes here
%   Temporary before creating a more generalizable one

arguments
    data
    spec_params
    base
    base_type = 1   % 1 = Troopa_Beh_Ephys_Hexa 220315 rn only
    to_plot logical = true
    x_tics_n = 8 % number of ticks on the x axis
end

% out
all_formattted = [];

% need to fix this so it's more generalizable with the data type
num_of_reaches = size(data,1);
out_all = cell(num_of_reaches,1);

% plotting
t_used = [-50, 200];

for i = 1:num_of_reaches  % per reach
    for j = 1:length(spec_params.fks)   % add more frequencies to test
        % neural data
        segment = sq(data(i,:,:));

        % base !! warning only for Troopa_hexa 220315 002 rn
        switch base_type
            case 1
                segment_base = sq(base(i,:,1:450)); % just before target comes on
            case 2
                disp('Add more LOL')
        end

        % initializing
        spec_out = cell(size(segment,1),1);
        raw_time = cell(size(segment,1),1);
        raw_freq = cell(size(segment,1),1);

        parfor iCh = 1:size(segment,1)
            % tfspec
            [spec_ch, f, ~] = tfspec(segment(iCh,:), ...
                                     [spec_params.tapers{j}, spec_params.smoothing{j}], ...
                                      spec_params.FS{j}, ...
                                      spec_params.w_step{j}, ...   
                                      spec_params.fks{j}, ...
                                      spec_params.pad{j},[],[], ...
                                      spec_params.flag{j});    % calc SPEC by pooling acrosss channels/trials
            % tfspec_base
            [spec_ch_base, f_base, ~] = tfspec(segment_base(iCh,:), ...
                                     [spec_params.tapers{j}, spec_params.smoothing{j}], ...
                                      spec_params.FS{j}, ...
                                      spec_params.w_step{j}, ...   
                                      spec_params.fks{j}, ...
                                      spec_params.pad{j},[],[], ...
                                      spec_params.flag{j});    % calc SPEC by pooling acrosss channels/trials
            
            % get calibration values
            [mu, sigma, ~, ~] = zlogcalib_alex_baselines(spec_ch_base, f_base);

            % z_log
            [zlog, ~] = zlogECoG_alex(spec_ch, f, mu, sigma);
            spec_out{iCh} = nanmean(zlog,2);
            
            % raw
            raw_time{iCh} = nanmean(log(spec_ch));  % across time
            raw_freq{iCh} = nanmean(log(spec_ch),2); % across freq

        end % parfor
        
        % formatting
        out = zeros(244,size(spec_out{1},1));
        out_raw = zeros(244,size(raw_time{1},2));
        out_raw2 = zeros(244, size(raw_freq{1},1));
        for l = 1:244
            out(l, :) = spec_out{l};
            out_raw(l, :) = raw_time{l};
            out_raw2(l, :) = raw_freq{l};
        end
        
        % save
        out_all{i} = out;

        % figures
        if to_plot
            n = x_tics_n; % # of ticks on the x-axis
            fig = figure(); set(fig, 'Visible', 'on'); set(gcf, 'Color', 'w');
            t = tiledlayout(3,1); title(t, sprintf('Reach %i', i))
            nexttile; 
            imagesc(out_raw);title('Raw HG - Log Mean (across time)');xlabel('Freq (Hz)');ylabel('Channel');colorbar;
            l_freq = get(gca, 'XLim');set(gca, 'XTick', linspace(l_freq(1), l_freq(2), 10)); set(gca, 'XTickLabel', linspace(spec_params.fks{1}(1), spec_params.fks{1}(2), 10))
            nexttile;
            imagesc(out_raw2);title('Raw HG - Log Mean (across freq)'); xlabel('Time (ms)');ylabel('Channel');colorbar;
            l_t = get(gca, 'XLim'); set(gca, 'XTick', linspace(l_t(1), l_t(2), n));set(gca, 'XTickLabel', linspace(t_used(1), t_used(2), n))
            nexttile;
            imagesc(out); xlabel('Time (ms)');ylabel('Channel');colorbar;
            l_t = get(gca, 'XLim'); set(gca, 'XTick', linspace(l_t(1), l_t(2), n));set(gca, 'XTickLabel', linspace(t_used(1), t_used(2), n))
            title(sprintf('HG Spec Z-log: w-step: %.3fsec taper: %.2fsec sm: %dHz', spec_params.w_step{j}, spec_params.tapers{j}, spec_params.smoothing{j}));
            input('Pres key for next reach')
        end
    end
end

all_formatted = zeros(num_of_reaches, size(out_all{1},1), size(out_all{1},2));
for l = 1:num_of_reaches
    all_formatted(l,:,:) = out_all{l};
end

end

