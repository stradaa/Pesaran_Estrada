function [fig, variance_values, mean_values] = plot_sigma_mu_overtime(data, chunk_size, name_of_plot)
%PLOT_SIGMA_MU_OVERTIME Plots mean and variance values over time
%
% Good for visualizing how the mean and variance change with time of the
% study.
%
% Alex 9-25-24

arguments
    data double % expects ch x time 2D double 
    chunk_size double
    name_of_plot string = "Mean and Variance"
end

% TODO reshape data if it is in 3D ch x time x trial or something

num_channels = size(data, 1);  % Number of channels (244)
num_timepoints = size(data, 2);  % Total number of time points (80800)
num_chunks = floor(num_timepoints / chunk_size);  % Number of chunks

% Initialize arrays to store mean and variance
mean_values = zeros(num_channels, num_chunks);
variance_values = zeros(num_channels, num_chunks);

% Loop over each channel and calculate the mean and variance for each chunk
for channel = 1:num_channels
    for chunk = 1:num_chunks
        % Define the indices for the current chunk
        start_idx = (chunk - 1) * chunk_size + 1;
        end_idx = chunk * chunk_size;
        
        % Extract the chunk of data for the current channel
        chunk_data = data(channel, start_idx:end_idx);
        
        % Compute the mean and variance of the chunk
        mean_values(channel, chunk) = mean(chunk_data);
        variance_values(channel, chunk) = var(chunk_data);
    end
end

% Generate time axis for plotting (one time point per chunk)
time_axis = (1:num_chunks) * chunk_size;

% Plot the mean values for each channel
fig = figure();t = tiledlayout(2,2);
fig.Position(3:4) = [1400 600];
nexttile(t,1);plot(time_axis, mean_values);title('Mean - line plot');xlabel('Trial #');ylabel('Mean');grid on;
nexttile(t,2);imagesc(mean_values);colorbar;title('Mean - imagesc');xlabel('Trial #');ylabel('Channel');
nexttile(t,3);plot(time_axis, variance_values');title('Variance - line plot');xlabel('Trial #');ylabel('Variance');
nexttile(t,4);imagesc(variance_values);title('Variance - imagesc');xlabel('Trial #');ylabel('Channel');colorbar;
title(t,name_of_plot);
set(fig, 'Color', 'w');

end

