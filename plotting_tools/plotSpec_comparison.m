function plotSpec_comparison(original_data, preprocessed_data, meta, show_plot, fs, freq_range, zoom_freq_range, window_length, overlap, corr_lags, title_text)
    %plotSpec_comparison
    %   Creates spectrogram and spectrum plots for a 2D input shape ch x t.
    %
    % Alex Estrada - Nov 15, 2024

    arguments
        original_data           double              % original_data: ch x t matrices
        preprocessed_data       double              % preprocessed_data: ch x t matrices
        meta        
        show_plot               logical = 1         % Visible = true by default
        fs                (1,1) double = 1000       % Sampling rate in Hz
        freq_range        (1,2) double = [0.1, 250] % Freq. range [min, max] (Hz)
        zoom_freq_range   (1,2) double = [0,250]    % Zoom range [min, max](Hz)
        window_length     (1,1) double = 0.1        % Length for spectrogram in (s)
        overlap           (1,1) double = 0.5        % 50% overlap
        corr_lags         (1,1) double = 30         % number of laggs for autocorrelation
        title_text              {string, char} = 'Comparison of Original and Preprocessed Data'
    end

    % Define the number of channels and time frames
    [ch, t] = size(original_data);
    time = (0:t-1) / fs;

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
    mean_spectrogram_orig = compute_mean_spectrogram(original_data);
    mean_spectrogram_prep = compute_mean_spectrogram(preprocessed_data);

    % Compute power spectra
    % compute_power_spectrum = @(data) ...
    %     deal(mean(abs(fft(data, [], 2)).^2, 1), ...
    %          std(abs(fft(data, [], 2)).^2, 0, 1) / sqrt(ch));

    spectrum_orig = cell(ch,1);
    spectrum_prep = cell(ch,1);
    for i = 1:ch
        [Pxx_1,freq_vector] = periodogram(original_data(i,:),rectwin(t),t,fs,"onesided", 'psd');
        [Pxx_2,~] = periodogram(preprocessed_data(i,:), rectwin(t),t,fs,"onesided", 'psd');

        spectrum_orig{i} = Pxx_1;
        spectrum_prep{i} = Pxx_2;
    end
    
    spectrum_orig = cell2mat(spectrum_orig);
    spectrum_prep = cell2mat(spectrum_prep);
    spectrum_orig = reshape(spectrum_orig, [],ch);
    spectrum_prep = reshape(spectrum_prep, [],ch);

    % mean power spectra
    mean_spectrum_orig = mean(spectrum_orig,2);
    mean_spectrum_prep = mean(spectrum_prep,2);

    % sd power spectra
    sd_spectrum_orig = std(spectrum_orig');
    sd_spectrum_prep = std(spectrum_prep');
        
    %% autocorrelation
    corr_orig = cell(ch,1);
    corr_prep = cell(ch,1);
    for i = 1:ch
        [acf_or, ~, ~] = autocorr(original_data(i,:),numLags=corr_lags);
        [acf_pr, lag_out, out_bounds] = autocorr(preprocessed_data(i,:),numLags=corr_lags);

        corr_orig{i} = acf_or;
        corr_prep{i}= acf_pr;
    end
    
    corr_orig = cell2mat(corr_orig);
    corr_prep = cell2mat(corr_prep);
    % mean corr
    mean_corr_orig = mean(corr_orig);
    mean_corr_prep = mean(corr_prep);
    % sd corr
    % sd_corr_orig = std(corr_orig);
    % sd_corr_prep = std(corr_prep);
    
    %% Plots
    % Initialize figure / plot
    fig = figure; tiledlayout(3, 2);

    %% Power spectrogram
    % Original
    nexttile;
    imagesc(time, frequencies, 10*log10(mean_spectrogram_orig));axis xy;
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('Original Power Spectrogram');colorbar;

    % Preprocessed
    nexttile;
    imagesc(time, frequencies, 10*log10(mean_spectrogram_prep));
    axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('Preprocessed Power Spectrogram');colorbar;
    l_freq = get(gca, 'XLim');
    n = ceil(diff(meta(1).bn)/100)+1;
    set(gca, 'XTick', linspace(l_freq(1), l_freq(2), n)); 
    set(gca, 'XTickLabel', linspace(meta(1).bn(1), meta(1).bn(2), n))
    
    %% power spectrum with SD
    nexttile;hold on;
    % Shaded area for standard deviation
    fill([freq_vector; flipud(freq_vector)], ...
         [mean_spectrum_orig - sd_spectrum_orig'; flipud(mean_spectrum_orig + sd_spectrum_orig')], ...
         'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shading
    % Plot mean spectrum
    plot(freq_vector, mean_spectrum_orig, 'b', 'LineWidth', 1.5); % Mean line
    hold off;xlabel('Frequency (Hz)');ylabel('Power');xlim(zoom_freq_range);
    title('Original Power Spectrum');grid on;
    
    % Plot preprocessed power spectrum with standard deviation shading
    nexttile;hold on;
    % Shaded area for standard deviation
    fill([freq_vector; flipud(freq_vector)], ...
         [mean_spectrum_prep - sd_spectrum_prep'; flipud(mean_spectrum_prep + sd_spectrum_prep')], ...
         'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shading
    % Plot mean spectrum
    plot(freq_vector, mean_spectrum_prep, 'r', 'LineWidth', 1.5); % Mean line
    hold off;xlabel('Frequency (Hz)');ylabel('Power');xlim(zoom_freq_range);
    title('Preprocessed Power Spectrum');grid on;

    %% ACF
    nexttile;
    stem(lag_out,mean_corr_orig);xlabel('Lag'); ylabel('\rho(k)');
    hold on; 
    h = line(lag_out,out_bounds(1)*ones(length(mean_corr_orig),1));
    h1 = line(lag_out,out_bounds(2)*ones(length(mean_corr_orig),1));
    set(h,'color',[1 0 0]);
    set(h1,'color',[1 0 0]);
    title('Autocorrelation - Original');
    % fill([lag_out'; flipud(lag_out)'], ...     % fill doesn't really work well
    %      [mean_corr_orig - sd_corr_orig'; flipud(mean_corr_orig + sd_corr_orig')], ...
    %      'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shading
    
    nexttile
    stem(lag_out,mean_corr_prep);xlabel('Lag'); ylabel('\rho(k)');
    hold on; 
    h = line(lag_out,out_bounds(1)*ones(length(mean_corr_prep),1));
    h1 = line(lag_out,out_bounds(2)*ones(length(mean_corr_prep),1));
    set(h,'color',[1 0 0]);
    set(h1,'color',[1 0 0]);
    title('Autocorrelation - Preprocessed');
    
    %% Final formatting
    sgtitle(title_text);
    set(fig, 'Color', 'w')

    if show_plot
        set(fig, 'Visible', 'on')
    end
end
