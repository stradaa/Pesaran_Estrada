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

globalMin = min(spec_avg(:));
globalMax = max(spec_avg(:));

% A small epsilon so the rightmost or topmost points aren’t out-of-range
xEdges(end) = xEdges(end) + 1e-9;
yEdges(end) = yEdges(end) + 1e-9;

ix = discretize(x, xEdges);  % which column (1..16)
iy = discretize(y, yEdges);  % which row    (1..16)

% Now create the tiled layout
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
end

colormap jet;

