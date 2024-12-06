function [S, F, T] = computeSpectrogram(signal, fs, smoothing, step_ms)
    % computeSpectrogram - Computes the spectrogram of a 1D signal.
    %
    %   spectrogram function in MATLAB computes the Short-Time Fourier 
    %   Transform (STFT) of the signal, which inherently produces 
    %   complex values because it represents both the amplitude and phase 
    %   information of the signal in the frequency domain. Here we only
    %   care about the power, so we compute the squared magnitude of the 
    %   given spectrogram.
    %
    % Inputs:
    %   signal     - Input signal (1xN double)
    %   fs         - Sampling frequency (Hz)
    %   smoothing  - Smoothing frequency (Hz)
    %   step_ms    - Window step size (ms)
    %
    % Outputs:
    %   S          - Spectrogram magnitude (frequency vs time)
    %   F          - Frequency bins (Hz)
    %   T          - Time bins (s)

    % Ensure the signal is a row vector
    if size(signal, 1) > 1
        signal = signal(:)';
    end
    
    % Define parameters
    win_length = fs / smoothing;           % Window length for smoothing
    step_size = round(step_ms / 1000 * fs); % Step size in samples
    overlap = win_length - step_size;     % Overlap between windows

    % Compute the spectrogram
    [S, F, T] = spectrogram(signal, win_length, overlap, [], fs);

    % Convert to power
    S = abs(S).^2; % Power of the spectrogram


    % Plot the spectrogram
    % figure;
    % imagesc(T, F, 10*log10(S));
    % axis xy;
    % xlabel('Time (s)');
    % ylabel('Frequency (Hz)');
    % title('Spectrogram');
    % colormap jet;
    % colorbar;
end
