function out = pos2matrix(data, electrodes, to_plot, w_labels, title_plot, color_lims)
%POS2MATRIX Proper electrode assignment
%   This function takes in 2D data and returns a correctly oriented movie
%   using the electrode positions. Returns a 3D squared array that can be
%   quickly animated with MovieSlider.m
%
%
% Alex Estrada - 06/18/2024- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data (244,:)    double          % expected size (channels x time)
    electrodes      struct          % electrodes (1x244 with 15 fields)
    to_plot        logical = true   % on by default
    w_labels       logical = true   % text with channel labels
    title_plot      string = ""     % name of the figure
    color_lims      double = []     % explicit color limits
end

%% Set-up
num_time_points = size(data, 2);
out = zeros(16,16, num_time_points); % initialize

%% General layout
x = zeros(244,0);
y = zeros(244,0);
elecRowCol = zeros(244,3);
for i = 1:244
    x(i) = electrodes(i).position.x;
    y(i) = electrodes(i).position.y;
    elecRowCol(i,2) = electrodes(i).position.row;
    elecRowCol(i,3) = electrodes(i).position.col;
    elecRowCol(i,1)= i;
end

%% Plot position
if to_plot
    fig = figure('Position', [671 441 1345 596]);
    set(fig, 'Color', 'w'); set(gcf, 'Visible', 'on')
    
    t = tiledlayout(1,2);
    if ~strcmp(t, "")
        title(t, title_plot)
    end
    nexttile
    scatter(x,y); axis padded;
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
    text(x,y,strsplit(num2str(1:244)))
    title('Position of Electrodes')
end

%% from agrita's getECoG_lmp_Layout.m ---
%fill the grid
maxRow=max(elecRowCol(:,2));

grid=zeros(16,16);
for ii=1:maxRow+1
    tmpRowsEle=find(maxRow-(ii-1)==elecRowCol(:,2)); % 
    
    for col=1:length(tmpRowsEle)
        tmpInd=tmpRowsEle(col);
        
        rowInd=elecRowCol(tmpInd,2)+1;
        colInd=elecRowCol(tmpInd,3)+1;
        elecNum=elecRowCol(tmpInd,1);
       
        grid(rowInd,colInd)=elecNum;
     end
end

gridLayout=flip(grid); % doing upside down - taking care of matlab and graph opps. convention

%% Assign and Reshape
% Iterate through the gridLayout
for row = 1:size(gridLayout, 1)
    for col = 1:size(gridLayout, 2)
        channel_index = gridLayout(row, col);
        if channel_index ~= 0
            out(row, col, :) = data(channel_index, :);
        end
    end
end

%% plot
if to_plot
    [numRows, numCols] = size(gridLayout);
    x_text = repmat(1:numCols,numCols,1); % generate x-coordinates
    y_text = repmat(1:numRows,numRows,1)'; % 
    
    textString = cellfun(@num2str,num2cell(gridLayout),'UniformOutput',false);
    
    nexttile
    imagesc(out(:,:,1));colorbar;
    if w_labels
        text(x_text(:), y_text(:), textString, 'HorizontalAlignment', 'Center');
    end
    title('Relative Tiled Layout')
    if ~isempty(color_lims)
        clim(color_lims)
    end
    set(gca, 'XTickLabel', [], 'YTickLabel', [])
end

end

