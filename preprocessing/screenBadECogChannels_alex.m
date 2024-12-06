function [out,L,U,C,sd] = screenBadECogChannels_alex(data,to_plot,RF,RM,data_y)
%SCREENBADECOGCHANNELS_ALEX Summary of this function goes here
%   Detailed explanation goes here

    arguments
        data double                 % ch x time
        to_plot logical = false     % generates the plot and makes it visible
        RF double = 1.5             % rejection factor
        RM {char, string} = 'rms'   % rejection method
        data_y {char, string} = 'Channel' % y-axis
    end

    % assuming data is already channels x time
    [a, b, c] = size(data);
    if c > 1
        warning('Attempting to reshape to channels x time (2D)')
        %reshape into channels x time
        data = reshape(permute(data,[2 1 3]), [b a*c]);
    end
    
    %% First Pass
    % want to look at the mean activity across time first
    mean_data = nanmean(data,2);
    [~,~,FPoutlier, L1,U1,C1] = rmoutliers(mean_data, 'mean', 'ThresholdFactor', 2);


    %% Main pass
    switch RM
        case 'rms' % first order rms
            temp = diff(data,2,2);
            sd = rms(temp,2);
        case 'mad'
        case 'sd'
            sd = nanstd(data,[],2);
        case 'kurtosis'
            % Quantifies prevalence of extreme values (outliers) in the 
            % data. Kurtosis is the fourth standardized moment of the data.
            % High values -> signal has large spikes or abrupt jumps.
            % Low values -> singal is more uniform

            temp = diff(data,2,2);
            sd = kurtosis(temp, [], 2); % across time
    end
    
    % get outlier
    [~,~,TFoutlier,L,U,C] = rmoutliers(sd, 'mean', 'ThresholdFactor', RF);
    
    %% Out (Aggregate)
    out = (FPoutlier | TFoutlier);
    
    % plotting
    if to_plot
        fig = figure; t = tiledlayout(2,1);title(t, 'Screening Overview');
        nexttile(1)
        plot(mean_data, '.'); hold on;
        plot(find(FPoutlier), mean_data(FPoutlier), 'r*')
        yline([L1 U1 C1], ":", ["Lower Threshold", "Upper Threshold", "Center Value"])
        xlabel(data_y); ylabel('Mean');
        title('First Pass | Mean value across time');
        hold off;

        nexttile(2)
        plot(sd, '.'); hold on;
        plot(find(TFoutlier), sd(TFoutlier), 'r*')
        yline([L U C], ":", ["Lower Threshold", "Upper Threshold", "Center Value"])
        ylabel(RM); xlabel(data_y);
        title(sprintf('%s Distribution | Method: %s | Factor: %.2f', data_y, RM, RF));
        hold off;

        set(fig, 'Color', 'w', 'Visible', 'on');
    end
    
    % output text
    id_bad = find(out);
    repstr = {repmat('%d ', 1, sum(out))};
    fprintf(['Bad Channels Numbers - screenBadECog:\n' repstr{:}, ...
             '\nTotal: %d\n'],id_bad , sum(out));

    test_thresh = 0;
    if test_thresh
        % manual selection here
        num_traces = 8; 
        thresh = 24;
        bn = [-600,1000];

        figure();tiledlayout(1,2);
        
        high_rms = find(sd>thresh);
        low_rsm = find(sd<thresh);
        
        shuffle_idx = randperm(length(high_rms));
        shuffle_idx2 = randperm(length(low_rsm));
        
        pick_idx = high_rms(shuffle_idx(1:num_traces));
        pick_idx2 = low_rsm(shuffle_idx2(1:num_traces));
        
        chs = data(pick_idx,:);
        chs2 = data(pick_idx2,:);
        
        nexttile(1);
        plot(chs');y1 = get(gca, 'YLim');
        legend(split(num2str(pick_idx(1:num_traces)')));
        title(['High ', RM]);
        l_freq = get(gca, 'XLim');
        n = ceil(diff(bn)/100)+1;
        set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
        set(gca, 'XTickLabel', linspace(bn(1), bn(2), n))
        
        nexttile(2);
        plot(chs2');ylim(y1)
        legend(split(num2str(pick_idx2(1:num_traces)')));
        title(['Low ', RM]);
        set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
        set(gca, 'XTickLabel', linspace(bn(1), bn(2), n))

        set(gcf, 'Color', 'w');
    end
end
