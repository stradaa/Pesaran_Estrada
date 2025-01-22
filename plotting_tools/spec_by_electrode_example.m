%% Loading data
load('C:\Users\Alex\Documents\Academics\Penn\PesaranLab\toy_data\electrodes.mat') % electrodes
load('C:\Users\Alex\Documents\Academics\Penn\PesaranLab\toy_data\spec_norm_data.mat') % spec

%% Average spec across trials
spec_avg = mean(spec_norm,4);
% now we have a 130 (time) x 128 (freq) x 244 (channel) average

%% Getting general electrode position
x = zeros(244,0); y = zeros(244,0);
elecRowCol = zeros(244,3);
for i = 1:244
    elecRowCol(i,1)= i;
    x(i) = electrodes(i).position.x; y(i) = electrodes(i).position.y;
    elecRowCol(i,2) = electrodes(i).position.row; elecRowCol(i,3) = electrodes(i).position.col;
end

figure();
scatter(x,y);axis padded;
set(gca, 'XTickLabel', [], 'YTickLabel', [])
text(x,y,strsplit(num2str(1:244)));
title('Electrode positions (X,Y)');

%%
% Suppose x, y are your 1xN doubles, N=244 in your case.

nCols = 16;
nRows = 16;

xEdges = linspace(min(x), max(x), nCols+1);
yEdges = linspace(min(y), max(y), nRows+1);

% A small epsilon so the rightmost or topmost points aren’t out-of-range
xEdges(end) = xEdges(end) + 1e-9;
yEdges(end) = yEdges(end) + 1e-9;

ix = discretize(x, xEdges);  % which column (1..16)
iy = discretize(y, yEdges);  % which row    (1..16)

% fix issues surrounding the holes
iy([30,27]) = 15;
iy([149,152]) = 2;

%% Now create the tiled layout
figure('Color','w');
t = tiledlayout(nRows, nCols, 'TileSpacing','none', 'Padding','none');
title(t, 'Discretized Tile-mapped Spectrograms');

% spec_avg is 130x128xN, one 130x128 image per electrode
for i = 1:numel(x)
    
    % row (1..16), but we want row=1 at top if we follow typical subplot order
    % If you want "largest y" at top, do row = nRows - iy(i) + 1
    row = nRows - iy(i) + 1;
    
    col = ix(i);  % 1..16
    
    % Convert to tile index (moving L->R, top->bottom)
    tile_index = (row - 1)*nCols + col;
    
    % If tile_index is NaN or out of range, skip or handle
    if isnan(tile_index) || tile_index<1 || tile_index>(nCols*nRows)
        continue
    end
    
    nexttile(tile_index);
    imagesc(spec_avg(:,:,i)');axis xy;

    clim([globalMin, globalMax]);

    axis off; axis tight;

    % if tile_index == 8
    %     disp('for testing')
    % end
end

colormap jet;

%%  Testing with a rotations of x90 degrees
globalMin = min(spec_avg(:));
globalMax = max(spec_avg(:));

% Now create the tiled layout
figure('Color','w', 'Position', [20 20 1000 1000]);
t = tiledlayout(nRows, nCols, 'TileSpacing','none', 'Padding','none');
title(t, 'Discretized Tile-mapped Spectrograms');

angle = 180;  % choose 0, 90, 180, or 270

for i = 1:numel(x)
    
    % 1) Convert (iy, ix) into "row" and "col" in your top-left coordinate system
    row = nRows - iy(i) + 1;  % top=1, bottom=16
    col = ix(i);              % left=1, right=16
    
    % 2) Rotate (row, col) by "angle" degrees
    [rowRot, colRot] = rotate_rc(row, col, nRows, nCols, angle);
    
    % 3) Compute tile index for the rotated row/col
    tile_index = (rowRot - 1)*nCols + colRot;
    
    % 4) If tile_index is out of range, skip. (Shouldn't happen if 1<=rowRot,colRot<=16)
    if tile_index < 1 || tile_index > nRows*nCols
        continue
    end
    
    nexttile(tile_index);
    imagesc(spec_avg(:,:,i)');
    axis xy; axis off; axis tight;
    clim([globalMin, globalMax]);
end

colormap jet;


%% Functions
function [rowNew, colNew] = rotate_rc(row, col, nRows, nCols, angle)
    % Rotate (row, col) on an nRows-by-nCols grid by "angle" degrees clockwise
    switch angle
        case 0
            rowNew = row;
            colNew = col;
        case 90
            % (r, c) -> (c, nRows-r+1)
            rowNew = col;
            colNew = nRows - row + 1;
        case 180
            % (r, c) -> (nRows-r+1, nCols-c+1)
            rowNew = nRows - row + 1;
            colNew = nCols - col + 1;
        case 270
            % (r, c) -> (nCols-c+1, r)
            rowNew = nCols - col + 1;
            colNew = row;
        otherwise
            error('Angle must be 0, 90, 180, or 270.');
    end
end

