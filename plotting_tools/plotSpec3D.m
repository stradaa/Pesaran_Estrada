function plotSpec3D(baseline_data, preprocessed_data, meta, show_plot, fs, freq_range, zoom_freq_range, window_length, overlap)
    %plotSpec_comparison
    %   Creates spectrogram and spectrum plots for a 3D input shape trial x ch x t.
    %
    % Alex Estrada - Nov 15, 2024

    arguments
        baseline_data           double              % baseline_data: trial x ch x t matrices
        preprocessed_data       double              % preprocessed_data: trial x ch x t matrices
        meta        
        show_plot               logical = 1         % Visible = true by default
        fs                (1,1) double = 1000       % Sampling rate in Hz
        freq_range        (1,2) double = [0.1, 250] % Freq. range [min, max] (Hz)
        zoom_freq_range   (1,2) double = [0,250]    % Zoom range [min, max](Hz)
        window_length     (1,1) double = 0.1        % Length for spectrogram in (s)
        overlap           (1,1) double = 0.25       % 25% overlap
    end

    % Define the number of channels and time frames
    [k, ch, t] = size(preprocessed_data);
    t_base = size(baseline_data,3);
    time = (0:t-1) / fs;
    time_base = (0:t_base-1)/fs;

    % Frequency vector for spectrogram
    frequencies = freq_range(1):1:freq_range(2);

    % Function to compute power spectrogram
    compute_mean_spectrogram = @(data) mean(cell2mat(arrayfun(@(ch) ...
        abs(spectrogram(data(ch, :), ...
        round(window_length*fs), ...
        round(overlap*window_length*fs), ...
        frequencies, fs, 'yaxis')).^2, ...
        1:ch, 'UniformOutput', false)), 3);

    % Compute mean power spectrograms
    spec_prep = cell(k,1);
    spec_base = cell(k,1);
    for i = 1:k
        temp_base = sq(baseline_data(i,:,:));
        temp_prep = sq(preprocessed_data(i,:,:));

        spec_base{i} = compute_mean_spectrogram(temp_base); % mean across ch
        spec_prep{i} = compute_mean_spectrogram(temp_prep); % mean across ch
    end
    
    % Reformat
    temp_base = cell2mat(spec_base');
    temp_prep = cell2mat(spec_prep');
    spec_base_out = reshape(temp_base, length(frequencies), length(spec_base{1}), k);
    spec_prep_out = reshape(temp_prep, length(frequencies), length(spec_prep{1}), k);
    
    % baseline
    [mu, sigma, ~, loga_freq] = zlogcalib_alex_baselines(temp_base', frequencies);
    [zlog, loga] = zlogECoG_alex(temp_prep', frequencies, mu, sigma);
    norm_prep = reshape(zlog', '')

    figure(); tiledlayout(2,2);
    nexttile();imagesc(log(sq(mean(spec_prep_out,3))));title('Log Mean Spec across ch and reach - Reach');
    axis xy;colorbar;xlabel('Time points');ylabel('Freq (Hz');
    nexttile();imagesc(log(sq(mean(spec_base_out,3))));title('Log Mean Spec across ch and reach - Baseline');
    axis xy;colorbar;xlabel('Time points');ylabel('Freq (Hz');
    nexttile();plot(loga_freq);title('Log mean baseline freq');xlabel('Freq (Hz)');ylabel('Log(dB)');
    
    nexttile();



    % norm_spec = cell2mat(spec_prep');
    % norm_spec = (norm_spec - mu)./sigma;
    % norm_spec = reshape(norm_spec, length(frequencies), length(spec_prep{1}), k);
    
    

    % Plotting the results (optional)
    if show_plot
        fig = figure; t = tiledlayout(1,2);

        temp_reach = 1; % temp
        % Prep
        nexttile(1,[1,2]);
        imagesc(10*log10(spec_prep{temp_reach}));
        axis xy;colorbar;title(sprintf('Window: %.2f, Overlap: %.2f', window_length, overlap))
        xlabel('Time (s)');ylabel('Frequency (Hz)');
        l_freq = get(gca, 'XLim');
        n = ceil(diff(meta(1).bn)/100)+1;
        set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
        set(gca, 'XTickLabel', linspace(meta(1).bn(1), meta(1).bn(2), n))

        % Preprocessed Data Spectrogram
        subplot(2, 1, 2);
        imagesc(time, frequencies, 10*log10(mean_spectrogram_prep));
        axis xy;colorbar;
        title(['Preprocessed Data: ' title_text]);
        xlabel('Time (s)');ylabel('Frequency (Hz)');
        l_freq = get(gca, 'XLim');
        n = ceil(diff(meta(1).bn)/100)+1;
        set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
        set(gca, 'XTickLabel', linspace(meta(1).bn(1), meta(1).bn(2), n))

        set(fig, 'Color', 'w')
        set(fig, 'Visible', 'on')
    end
end
