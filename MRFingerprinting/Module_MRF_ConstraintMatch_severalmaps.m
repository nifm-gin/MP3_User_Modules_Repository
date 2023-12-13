function [files_in,files_out,opt] = Module_MRF_ConstraintMatch_severalmaps(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

	%   % define every option needed to run this module
	% --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'dictionary_folder_filename',  'Dictionary Folder'};
    module_option(:,2)   = {'prefix',           'match_constraint_'};
    module_option(:,3)   = {'method',           'DBM'};
    module_option(:,4)   = {'filtered',         'No'};
    module_option(:,5)   = {'removed',          0};
    module_option(:,6)   = {'augment',          60};
    
    module_option(:,7)   = {'RefInput',         1};
    module_option(:,8)   = {'InputToReshape',   1};
    module_option(:,9)   = {'Table_in',         table()};
    module_option(:,10)  = {'Table_out',        table()};
    module_option(:,11)  = {'folder',           table()};
    module_option(:,12)  = {'OutputSequenceName','AllName'};
    module_option(:,13)  = {'clipValues',       ''};
    module_option(:,14)  = {'confMaps',         'No'};
    module_option(:,15)  = {'Params',           'Vf'};
    module_option(:,16)  = {'K',                50};
    module_option(:,17)  = {'Lw',               0};
    module_option(:,18)  = {'cstrS',            'd*'};
    module_option(:,19)  = {'cstrG',            ''};
    module_option(:,20)  = {'mkScoreMap', 'Yes'};
    module_option(:,21)  = {'mkMatchedVolume', 'Yes'};
    module_option(:,22)  = {'echoes', '1'};
    module_option(:,23)   = {'RestType', 'B1rel'};
    module_option(:,25)   = {'otherRest', 'No'};
    module_option(:,26)   = {'RestType2', 'Ttwo_star'};

    module_option(:,24)   = {'finalNorm', 'Yes'};

    
    opt.Module_settings  = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    
    %% list of everything displayed to the user associated to their 'type'
    % --> user_parameter(1,:) = user_parameter_list
    % --> user_parameter(2,:) = user_parameter_type
    % --> user_parameter(3,:) = parameter_default
    % --> user_parameter(4,:) = psom_parameter_list
    % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
    % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional.
    % --> user_parameter(7,:) = Help : text data which describe the parameter (it
    % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'The dictionary-based matching (DBM) reconstruction method is based on the paper: Ma, Dan, et al. "Magnetic resonance fingerprinting." Nature (2013)'
        'The dictionary-based learning (DBL) reconstruction method is based on the paper: Boux, Fabien, et al. [work in progress]'
        ''
        'Prerequisite:'
        '      - Put your ''PRE_*.json'' and ''POST_*.json'' dictionary files (MGEFIDSE pre and post simulated scans) in the ''MP3/data/dictionaries'' folder'
        '      - (or) Put your ''DICO.mat'' dictionary file ratio between the MGEFIDSE pre and post simulated scans in the ''MP3/data/dictionaries'' folder'
        '      - (or - DBL method only) Put your ''MODEL_*.mat'' model file'
        ''
        'The dictionaries are designed with the Mrvox simulation tool, see [.]'
        }'};
    
    user_parameter(:,2)   = {'Select the acquired data','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    folder_files	= dir(fullfile(s{1}, 'data/dictionaries/'));
    folder_files    = folder_files([folder_files.isdir]);
    opt.Module_settings.folder = fullfile(s{1}, 'data/dictionaries/');
    if isempty(folder_files), folder_files(1).name = ' '; end
    
    user_parameter(:,3)   = {'   .Dictionary Pre/Post file folder','cell', {folder_files.name}, 'dictionary_folder_filename','','Mandatory',...
        {'Select the folder containing Pre/Post dico files (.json), ratio dico file (.mat) and/or model file (.mat)'}};
    
    user_parameter(:,4)   = {'   .Prefix','char', '', 'prefix', '', '',...
        {'Choose a prefix for your output maps'}};
   

    user_parameter(:,8)   = {'   .Echoes to keep?','char','', 'echoes', '', '',...
        {'How to take the echoes during matching (all echoes, or only 1st echoes, etc...). It allows also to choose starting and ending point : start:step:end'}};
    
    user_parameter(:,5)   = {'   .Parameters','check', ...
        {'T2', 'B0off', 'T1', 'df', 'B1rel',  'Ttwo_star', 'gamma', 'Vf', 'SO2', 'R', 'VSI', 'DH2O', 'B0theta', 'khi', 'Hct', 'Label'},...
        'Params', '', '',...
        {'Select the parameters considered in the model'
        }'};
    user_parameter(:,6)   = {'   .Clip','char', '0 inf' 'clipValues', '', '',...
        {'Clip produced maps to desired values. Indicate "min(param 1), max(param 1); min(param 2), max(param 2); ...", default is no clip'}};
    user_parameter(:,9)   = {'   .Confidence Maps','cell', {'Yes','No'}, 'confMaps', '', '',...
        {'Produce confidence maps, default = No'}};
    
    user_parameter(:,10)   = {'   .Make ScoreMap?','cell', {'Yes','No'}, 'mkScoreMap', '', '',...
        {'Select ''Yes'' add score map in outputs  (recommanded ''Yes'')';''}'};
    user_parameter(:,7)   = {'   .produce Matched Volume?','cell', {'Yes','No'}, 'mkMatchedVolume', '', '',...
        {'Select ''Yes'' will produce a volume composed of the matched signals from the dico, L2 normalized (recommanded ''Yes'')'; 'WARNING: does not work with DBL'}'};
 
    user_parameter(:,11)   = {'Indicate the type of restriction','cell', {'B1rel', 'T1', 'T2', 'df', 'Ttwo_star'}, 'RestType', '', '',...
        {}'};
    user_parameter(:,12)   = {'Select the scans to use for restriction','XScan','','',{'SequenceName'}, '','A T1/T2/T2star map should be in ms, B1map normalized, df map in Hz'};
    user_parameter(:,13)   = {'Use a second restriction map?','cell', {'Yes','No'}, 'otherRest', '', '',...
        {'Select ''Yes'' will use two parameters restriction. '}'};
 
    user_parameter(:,14)   = {'   .IF YES : Indicate the type of second restriction','cell', {'B1rel', 'T1', 'T2', 'df', 'Ttwo_star'}, 'RestType2', '', '',...
        {}'};
    user_parameter(:,15)   = {'   .IF YES : Select the scans to use for second restriction','XScan','','',{'SequenceName'}, '','A T1/T2/T2star map should be in ms, B1map normalized, df map in Hz'};
    
    user_parameter(:,16)   = {'Normalize final signal', 'cell', {'Yes', 'No'}, 'finalNorm', '','Optional', 'Do you want to normalize (N2) the final combination before MRF?'};

    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)', 'VariableNames', VariableNames);
    %%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1    = {''};
    files_out.In1   = {''};
    return
    
    
end
%%%%%%%%


if isempty(files_out)
    
    opt.Params = opt.Params(cell2mat(opt.Params(:,2)),1);
    
    opt.Table_out = opt.Table_in(1,:);
    
    for i = 1:numel(opt.Params)
        opt.Table_out(i,:) = opt.Table_out(1,:);
        opt.Table_out(i,:).Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            opt.Table_out.SequenceName(i) = categorical(cellstr([char(opt.prefix), char(opt.Params{i})]));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            opt.Table_out.SequenceName(i) = categorical(cellstr([char(opt.Table_out.SequenceName(i)), opt.Params{i}]));
        end
        opt.Table_out.Filename(i) = categorical(cellstr([char(opt.Table_out.Patient(i)), '_', char(opt.Table_out.Tp(i)), '_', char(opt.Table_out.SequenceName(i))]));
        opt.Table_out.IsRaw(i) = categorical(cellstr('0'));
        files_out.In1{i} = [char(opt.Table_out.Path(i)), char(opt.Table_out.Filename(i)) '.nii'];
    end

    if strcmp(opt.mkScoreMap, 'Yes')
        scoreRow.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            scoreRow.SequenceName = categorical(cellstr([char(opt.prefix), 'ScoreMap']));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            scoreRow.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), 'ScoreMap']));
        end
        scoreRow.Filename = categorical(cellstr([char(opt.Table_out.Patient(1)), '_', char(opt.Table_out.Tp(1)), '_', char(scoreRow.SequenceName)]));
        scoreRow.IsRaw = categorical(cellstr('0'));
        
        scoreRow.Group = opt.Table_out.Group(1);
        scoreRow.Patient = opt.Table_out.Patient(1);
        scoreRow.Tp = opt.Table_out.Tp(1);
        scoreRow.Type = opt.Table_out.Type(1);
        opt.Table_out = [opt.Table_out; struct2table(scoreRow)];
        
        files_out.In1{end+1} = [char(opt.Table_out.Path(end)), char(opt.Table_out.Filename(end)) '.nii'];
    end
    
    if strcmp(opt.mkMatchedVolume, 'Yes')
        matchedVolRow.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            matchedVolRow.SequenceName = categorical(cellstr([char(opt.prefix), 'MatchedVol']));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            matchedVolRow.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), 'MatchedVol']));
        end
        matchedVolRow.Filename = categorical(cellstr([char(opt.Table_out.Patient(1)), '_', char(opt.Table_out.Tp(1)), '_', char(matchedVolRow.SequenceName)]));
        matchedVolRow.IsRaw = categorical(cellstr('0'));
        
        matchedVolRow.Group = opt.Table_out.Group(1);
        matchedVolRow.Patient = opt.Table_out.Patient(1);
        matchedVolRow.Tp = opt.Table_out.Tp(1);
        matchedVolRow.Type = opt.Table_out.Type(1);
        opt.Table_out = [opt.Table_out; struct2table(matchedVolRow)];
        
        files_out.In1{end+1} = [char(opt.Table_out.Path(end)), char(opt.Table_out.Filename(end)) '.nii'];
    end
    
end

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Smoothing:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.finalNorm = strcmp(opt.finalNorm, 'Yes');

if ~isempty(opt.clipValues)
    clipValues = str2num(opt.clipValues);
    assert(numel(opt.Params) == size(clipValues, 1), 'The number of clipping couples (min, max) provided does not match the number of parameters')
    flagClip = 1;
else
    flagClip = 0;
end
% Read json files and create ratio dictionary
opt.dictionary_folder_filename = fullfile(opt.folder, opt.dictionary_folder_filename);
d = dir(opt.dictionary_folder_filename);

n = d(contains({d.name}, 'DICO'));
opt.dictionary_filename = n.name;
% d   = dir(opt.dictionary_folder_filename);
% opt.dictionary_pre_filename     = d(contains({d.name}, 'PRE_'));
% if ~isempty(opt.dictionary_pre_filename)
%     opt.dictionary_pre_filename     = opt.dictionary_pre_filename.name;
% end
% opt.dictionary_post_filename    = d(contains({d.name}, 'POST_'));
% if ~isempty(opt.dictionary_post_filename)
%     opt.dictionary_post_filename    = opt.dictionary_post_filename.name;
% end


% Load and if necessary, smooth observations
% Xpost   = niftiread(files_in.In2{1});
% Xpre    = niftiread(files_in.In1{1});
scan_of_reference.header = spm_vol(files_in.In1{1});
scan_of_reference.data=  read_volume(scan_of_reference.header, scan_of_reference.header, 0, 'Axial');
Xobs = scan_of_reference.data;

json_filename   = split(files_in.In1{1},'.');
json_filename{2} = '.json';
Obs  = ReadJson([json_filename{1} json_filename{2}]);




% Generate conc signals from scans (and TODO: ROI if given)
% Xpre = Xpre ./ vecnorm(Xpre,2,4);
% Xpost = Xpost ./ vecnorm(Xpost,2,4);
% Xobs = cat(4, Xpre, Xpost);
Xobs = Xobs ./ vecnorm(Xobs,2,4);

% Dictionary filename
dico_filename = [opt.dictionary_folder_filename filesep 'DICO.mat'];

% If, dico exists, load it, else, create it and save it
if exist(dico_filename,'file')
    load(dico_filename)

else
%         Pre     = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
%         Post    = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_post_filename]);
D = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
Dico.MRSignals{1}       = abs(Pre.MRSignals);
Dico.MRSignals{2}       = abs(Post.MRSignals);
%         Dico.MRSignals      = abs(Post.MRSignals) ./ abs(Pre.MRSignals); % Ratio post/pre signals 
Dico.Tacq           = Pre.Sequence.Tacq;

Dico.Parameters.Par = Pre.Parameters.Par; % Parameters used to simulate X signals
Dico.Parameters.Labels = Pre.Parameters.Labels;
clear Pre Post
save(dico_filename,'Dico')
end

%% Select wanted echoes in the sequence
echoes_to_keep = str2num(opt.echoes);

% If necessary, choose echoes
if size(echoes_to_keep, 2) ~= 1
% crop the acquisition
    Xobs = Xobs(:,:,:, echoes_to_keep);
%     Obs.EchoTime.value = Obs.EchoTime.value(echoes_to_keep);
% crop the dico
    Dico.MRSignals{1,1} = Dico.MRSignals{1,1}(:, echoes_to_keep);
    Dico.Tacq = Dico.Tacq(echoes_to_keep);
%     clear tmp
end

% Permute
Xobs(Xobs == inf) = nan;
Xobs	= permute(Xobs, [1 2 4 3]);


%% DBM

%     if size(Xobs,length(size(Xobs)))/2 ~= size(Dico.MRSignals{1},2) % Divided by two because Xobs = [Xpre Xpost]
%         warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
%         Dico.MRSignals{1} = interp1(Dico.Tacq(1:size(Dico.MRSignals{1},2)), Dico.MRSignals{1}', Obs.EchoTime.value'*1e-3)';
%         Dico.MRSignals{2} = interp1(Dico.Tacq(1:size(Dico.MRSignals{2},2)), Dico.MRSignals{2}', Obs.EchoTime.value'*1e-3)';       
%     end
Dico.Tacq   = Obs.EchoTime.value'*1e-3;
%remove row containning nan values
if ~iscell(Dico.MRSignals)
    tmpSave = Dico.MRSignals;
    Dico.MRSignals = {};
    Dico.MRSignals{1} = tmpSave;
end
nn = ~any(isnan(Dico.MRSignals{1}),2);
Dico.MRSignals{1}  = Dico.MRSignals{1}(nn,:);
%     Dico.MRSignals{2}  = Dico.MRSignals{2}(nn,:);
Dico.Parameters.Par = Dico.Parameters.Par(nn,:);

%% Give dico its shape
% This step is already performed when interpolating the signal with
    % Obs.EchoTime.value, which is itself already cropped
% if opt.removed > 0
%     Dico.MRSignals{1}   = Dico.MRSignals{1}(:, 1 : end - opt.removed);
%     Dico.MRSignals{2}   = Dico.MRSignals{2}(:, 1 : end - opt.removed);
%     Dico.Tacq           = Dico.Tacq(1 : end - opt.removed);
% end

Dico.MRSignals{1}       = Dico.MRSignals{1} ./ vecnorm(Dico.MRSignals{1}, 2, 2);



%% Dico restriction
% Check input restriction maps
switch opt.RestType
    case 'T1'
        inPar = 'T1';
        Fact = 1e3;
    case 'T2'
        inPar = 'T2';
        Fact = 1e3;
    case 'Ttwo_star'
        inPar = 'T2';
        Fact = 1e3;
    case 'B1rel'
        inPar = 'B1rel';
        Fact = 1;
    case 'df'
        inPar = 'df';
        Fact = 1;
end

if strcmp(opt.otherRest, 'Yes')
    switch opt.RestType2
        case 'T1'
            inPar2 = 'T1';
            Fact2 = 1e3;
        case 'T2'
            inPar2 = 'T2';
            Fact2 = 1e3;
        case 'Ttwo_star'
            inPar2 = 'T2';
            Fact2 = 1e3;
        case 'B1rel'
            inPar2 = 'B1rel';
            Fact2 = 1;
        case 'df'
            inPar2 = 'df';
            Fact2 = 1;
    end
end

%% Reduce the number of search
colNb               = find(strcmp(Dico.Parameters.Labels,inPar));
colNb2               = find(strcmp(Dico.Parameters.Labels,inPar2));

RestMap.header =  spm_vol(files_in.In2{1});
RestMap.data = read_volume(RestMap(1).header, scan_of_reference.header, 0, 'Axial');
RestMap = RestMap.data;
Values = unique(Dico.Parameters.Par(:, colNb));
% Loop through each voxel in RestMap
for i = 1:numel(RestMap)
    % Find the nearest value in Values vector
    [~, idx] = min(abs(RestMap(i) - Values));    
    % Round the voxel to the nearest value in Values vector
    RestMap(i) = Values(idx);
end
if nnz(isnan(Values))
    warning('%i voxels will not be evaluated as their %s value falls outside the dictionary range', nnz(isnan(Values)), inPar)
end


if strcmp(opt.otherRest, 'Yes')
    RestMap2.header =  spm_vol(files_in.In3{1});
    RestMap2.data = read_volume(RestMap2(1).header, scan_of_reference.header, 0, 'Axial');
    RestMap2 = RestMap2.data;
    Values2 = unique(Dico.Parameters.Par(:, colNb2));
    % Loop through each voxel in RestMap
    for i = 1:numel(RestMap2)
        % Find the nearest value in Values vector
        [~, idx] = min(abs(RestMap2(i) - Values2));    
        % Round the voxel to the nearest value in Values vector
        RestMap2(i) = Values2(idx);
    end
    if nnz(isnan(Values2))
        warning('%i voxels will not be evaluated as their %s value falls outside the dictionary range', nnz(isnan(Values2)), inPar)
    end
end
% minDicoVal          = min(Dico.Parameters.Par(:,colNb))*Fact; % min of the dico
% maxDicoVal          = max(Dico.Parameters.Par(:,colNb))*Fact; % max of the dico
% Values(Values*(1+opt.RelErr) < minDicoVal) = nan;
% Values(Values*(1-opt.RelErr) > maxDicoVal) = nan;



% Dico Restriction
for v=1:numel(Values)
    for w=1:numel(Values2)
        % Removing dico entries where parameter of interest is out of
        % the tolerated range
        if isnan(Values(v)) %If value at this iteration is nan, don't consider it
            continue
        end
        % if not nan, get coordinates of corresponding voxels
        [row, col, sl] = ind2sub(size(RestMap), find(RestMap == single(Values(v)) & RestMap2 == single(Values2(w)))); % Get coordinates of voxels considered at this iteration

        % remove dico entries out of tolerated range
        toRemoveInf     = Dico.Parameters.Par(:,colNb)*Fact < Values(v) & Dico.Parameters.Par(:,colNb2)*Fact2 < Values2(w);
        toRemoveSup     = Dico.Parameters.Par(:,colNb)*Fact > Values(v) & Dico.Parameters.Par(:,colNb2)*Fact2 > Values2(w);
        toKeep          = ~(toRemoveInf + toRemoveSup);

        if nnz(toKeep) == 0
            warning('Value %i could not be evaluated as restricted dico is empty, %i voxels concerned\n', Values(v), numel(row))
            continue
        end

        % Copying the restricted dico
        TmpDico{1}.MRSignals = Dico.MRSignals{1,1}(toKeep, :);
        TmpDico{1}.Parameters.Par = Dico.Parameters.Par(toKeep, :);
        TmpDico{1}.Parameters.Labels = Dico.Parameters.Labels;

        localXobs = zeros(numel(row), size(Xobs, 3));
        for k=1:numel(row)
            localXobs(k,:) = Xobs(row(k), col(k), :, sl(k)); % Keep only the corresponding voxels in the observation
        end
        % TODO: find something nicer than this permute trick
        Estimation  = AnalyzeMRImages(localXobs,TmpDico,opt.method,[]);
        Map.Y       = permute(Estimation.GridSearch.Y, [1 2 4 3]);



    %%
        % Extract maps (and modify unit if necessary)
        count = 1;
        for i = 1:length(TmpDico{1}.Parameters.Labels)
            tmp = split(TmpDico{1}.Parameters.Labels{i},'.',2);
            if any(strcmp(tmp{end}, opt.Params))
                Labels{count} = tmp{end};      
                switch tmp{end}
                    %If the ression method is performed, extract also the confidence maps
                    case {'Vf', 'SO2'} %convert to percent
                        for j = 1:numel(row)
                            MapStruct{count}(row(j), col(j), sl(j))     = 100*Map.Y(j,:,:,i);
                        end
                    case {'VSI', 'R'} % convert m to µm
                        for j = 1:numel(row)
                            MapStruct{count}(row(j), col(j), sl(j))     = 1e6*Map.Y(j,:,:,i);
                        end
                    case {'T1', 'Ttwo_star', 'T2'} % convert s to ms
                        for j = 1:numel(row)
                            MapStruct{count}(row(j), col(j), sl(j))     = Map.Y(j,:,:,i);
                        end
                    otherwise
                        for j = 1:numel(row)
                            MapStruct{count}(row(j), col(j), sl(j)) 	= Map.Y(j,:,:,i);
                        end
                end   % end switch tmp
                count   = count +1;
            end
        end
        for j = 1:numel(row)
            ScoreMap(row(j), col(j), sl(j)) = Estimation.scoreMap(j);
            MatchedVol(row(j), col(j), sl(j), :) = Estimation.matchedVol(j,:);
        end
    end
end

if strcmp(opt.mkScoreMap, 'Yes')
    MapStruct{end+1} = ScoreMap;
    Labels{end+1} = 'ScoreMap';
end

if strcmp(opt.mkMatchedVolume, 'Yes')
    MapStruct{end+1} = MatchedVol; % direct map of the signals from dico 
    Labels{end+1} = 'MatchedVol';
end
    
% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J       = ReadJson(jsonfile);
J       = KeepModuleHistory(J, struct('files_in',files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 


% Reoriented and save nifti maps
% nifti_header = spm_vol(files_in.In1{1});
info = niftiinfo(files_in.In1{1});
info.Filemoddate = char(datetime('now'));

for i = 1:length(files_out.In1)
    
    for j = 1:numel(Labels)
        [~, fOut, ~] = fileparts(files_out.In1{i});  
        param_tmp = split(fOut, '_');
        if contains(param_tmp{end},Labels{j})
            [path, name, ~] = fileparts(files_out.In1{i});
            WriteJson(J, [path, '/', name, '.json'])
            info2 = info;
            info2.Filename = files_out.In1{i};
            info2.ImageSize = size(MapStruct{j});
            info2.PixelDimensions = info2.PixelDimensions(1:length(size(MapStruct{j})));
            info2.Datatype = class(MapStruct{j});
            MapStruct{j} = write_volume(MapStruct{j}, scan_of_reference.header, 'Axial');
            niftiwrite(MapStruct{j}, files_out.In1{i}, info2);
        end
    end
end