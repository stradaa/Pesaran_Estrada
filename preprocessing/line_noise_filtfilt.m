function filtered_data = line_noise_filtfilt(data, fs, bad_ch, show_plot)
%LINE_NOISE_FILTFILT Summary of this function goes here
%   Detailed explanation goes here

arguments
    data    % trials x channels x time
    fs      % sampling rate
    bad_ch  % channels to ignore
    show_plot
end

% Define parameters
harmonics = 60:60:round(fs/2); % Line noise harmonics up to Nyquist frequency (500 Hz)

% Design notch filters
filters = cell(length(harmonics), 1);
for i = 1:length(harmonics)
    filters{i} = designfilt('bandstopiir', 'FilterOrder', 2, ...
                            'HalfPowerFrequency1', harmonics(i)-1, ...
                            'HalfPowerFrequency2', harmonics(i)+1, ...
                            'SampleRate', fs);
end

filtered_data = data; % Initialize filtered data

% Apply filters sequentially to each trial and channel
for trial = 1:size(data, 1)
    for channel = 1:size(data, 2)
        if ~ismember(channel, bad_ch)
            signal = squeeze(data(trial, channel, :));
            for i = 1:length(filters)
                signal = filtfilt(filters{i}, signal);
            end
            filtered_data(trial, channel, :) = signal;
        end
    end
end

if show_plot
    % Select a specific channel and trial for plotting
    channel_index = 150;
    trial_index = 1;
    original_signal = squeeze(data(trial_index, channel_index, :));
    filtered_signal = squeeze(filtered_data(trial_index, channel_index, :));
    
    % Plot the spectrum before and after filtering
    nfft = 2^nextpow2(length(original_signal));
    frequencies = (0:(nfft/2)) * (fs / nfft);
    original_spectrum = abs(fft(original_signal, nfft)) / length(original_signal);
    filtered_spectrum = abs(fft(filtered_signal, nfft)) / length(filtered_signal);
    
    figure;
    subplot(2, 1, 1);
    plot(frequencies, log(2*original_spectrum(1:nfft/2+1)));
    title('Original Signal Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    
    subplot(2, 1, 2);
    plot(frequencies, log(2*filtered_spectrum(1:nfft/2+1)));
    title('Filtered Signal Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');

    set(gcf, 'Color', 'w', 'Visible', 'on')
end

end

