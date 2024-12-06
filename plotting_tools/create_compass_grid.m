function compass_grid = create_compass_grid(grid_size, to_plot)
% Create a grid and insert patterns that resemble the 
% shapes of "N", "S", "W", and "E".
%
% This function is used to determine the proper transformation for when
% getting the regristration all figured out for any plot.
%
% Alex 09-23-24

arguments
    grid_size = 150
    to_plot = false
end

compass_grid = zeros(grid_size, grid_size);

% Create binary masks for each letter (N, S, W, E)
% Binary pattern for 'N'
letter_N = [
    1 0 0 0 1;
    1 1 0 0 1;
    1 0 1 0 1;
    1 0 0 1 1;
    1 0 0 0 1;
    ];

% Binary pattern for 'S'
letter_S = [
    1 1 1;
    1 0 0;
    1 1 1;
    0 0 1;
    1 1 1;
    ];

% Binary pattern for 'W'
letter_W = [
    1 0 0 0 1;
    1 0 0 0 1;
    1 0 1 0 1;
    1 0 1 0 1;
    1 1 0 1 1;
    ];

% Binary pattern for 'E'
letter_E = [
    1 1 1;
    1 0 0;
    1 1 1;
    1 0 0;
    1 1 1;
    ];

% Scale the binary patterns to fit the grid (5x scaling factor)
scale_factor = 5;
letter_N = imresize(letter_N, scale_factor, 'nearest');
letter_S = imresize(letter_S, scale_factor, 'nearest');
letter_W = imresize(letter_W, scale_factor, 'nearest');
letter_E = imresize(letter_E, scale_factor, 'nearest');

% Get the size of the scaled letters
[n_height, n_width] = size(letter_N);
[s_height, s_width] = size(letter_S);
[w_height, w_width] = size(letter_W);
[e_height, e_width] = size(letter_E);

% Place the letters into the compass grid
% North ('N') on the top-center
compass_grid(1:n_height, round(grid_size/2)-round(n_width/2)+1:round(grid_size/2)-round(n_width/2)+n_width) = letter_N;

% South ('S') on the bottom-center
compass_grid(end-s_height+1:end, round(grid_size/2)-round(s_width/2)+1:round(grid_size/2)-round(s_width/2)+s_width) = letter_S;

% West ('W') on the center-left
compass_grid(round(grid_size/2)-round(w_height/2)+1:round(grid_size/2)-round(w_height/2)+w_height, 1:w_width) = letter_W;

% East ('E') on the center-right
compass_grid(round(grid_size/2)-round(e_height/2)+1:round(grid_size/2)-round(e_height/2)+e_height, end-e_width+1:end) = letter_E;

% Plot
if to_plot
    figure();
    imagesc(compass_grid);
    colormap(gray);
    axis off;
    axis image;
end

end

