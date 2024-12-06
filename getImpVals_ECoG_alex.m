function [electrodes, elec_imps, gridLayout, impGridLayout] = getImpVals_ECoG_alex(meta)
%GETIMPVALS_ECOG Creates plot of impedence values on ECoG given experiment
%   Adapted from:
%       [elec_imps]=getImpVals_ECoG_rhd(MONKEYDIR,day,experiment) (Agrita)
%
% Alex Estrada - 11/14/2024

arguments
    meta 
end

electrodes = [];
elec_imps = [];
gridLayout = [];
impGridLayout = [];

% check if multiple recs
if size(meta,2) > 1
    m = meta(1);        % defaults to first rec idx
else
    m = meta;
end

% prompt user for file (methods/naming of file are not standard yet)
[imp_file, path_name] = uigetfile([m.MONKEYDIR, '/', m.day, '/*.csv'], ...
                                'Get Impedance File');
if path_name == 0 
    disp('No electrode impedence CSV selected.')
    return
end
% get electrode imps
imp_files = {[path_name, imp_file]};
[elec_imps, ~] = match_impedances_from_intan_csv_to_elecs(imp_files, m.experiment, m.imd,  'intan_rhd');

% electrodes
electrodes = m.experiment.hardware.microdrive(1).electrodes;


%% adapted from getECoG_Imp_Layout() Trropa_Be_Ep_Hexa/matlab/mfiles/agrita
expDef_electrodes=electrodes;
numElecs = length(expDef_electrodes);
elecXY = zeros(numElecs,3);
elecRowCol = zeros(numElecs,3);

for i=1:numElecs
    tmp=expDef_electrodes(i).position;
    tmpElec=expDef_electrodes(i).acquisitionid;  % electrode label --> sensor
%   tmpElec=expDef_electrodes(i).channelid;      % channel  label --> daq

    elecXY(i,1)=tmpElec;
    elecXY(i,2)=tmp.x;
    elecXY(i,3)=tmp.y;
    
    elecRowCol(i,1)=tmpElec;
    elecRowCol(i,2)=tmp.row;
    elecRowCol(i,3)=tmp.col;
end

%fill the grid
maxRow=max(elecRowCol(:,2));

grid=zeros(16,16);
imp_grid=nan(16,16);
imp_grid_kohm=nan(16,16);

for ii=1:maxRow+1
    tmpRowsEle=find(maxRow-(ii-1)==elecRowCol(:,2)); % 
    
    for col=1:length(tmpRowsEle)
        tmpInd=tmpRowsEle(col);
        
        rowInd=elecRowCol(tmpInd,2)+1;
        colInd=elecRowCol(tmpInd,3)+1;
        elecNum=elecRowCol(tmpInd,1);
       
        grid(rowInd,colInd)=elecNum;
        imp_grid(rowInd,colInd)=elec_imps(elecNum); %ohms
        imp_grid_kohm(rowInd,colInd)= elec_imps(elecNum)./1e3; %kiloOhm
    end
end

gridLayout=flip(grid); % doing upside down - taking care of matlab and graph opps. convention
impGridLayout=flip(imp_grid_kohm); %kilo Ohms
end

