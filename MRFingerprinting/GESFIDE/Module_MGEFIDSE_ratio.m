function [files_in,files_out,opt] = Module_MGEFIDSE_ratio(files_in,files_out,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

	%   % define every option needed to run this module
	% --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'dictionary_folder_filename',  'Dictionary Folder'};
    module_option(:,2)   = {'prefix',           'DBL_'};
    module_option(:,3)   = {'method',           'DBL'};
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
    
    user_parameter(:,2)   = {'Select the MGEFIDSE Pre scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the MGEFIDSE Post scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    
    s               = split(mfilename('fullpath'),'MP3_pipeline',1);
    folder_files	= dir(fullfile(s{1}, 'data/dictionaries/'));
    folder_files    = folder_files([folder_files.isdir]);
    opt.Module_settings.folder = fullfile(s{1}, 'data/dictionaries/');
    if isempty(folder_files), folder_files(1).name = ' '; end
    user_parameter(:,4)   = {'   .Dictionary Pre/Post file folder','cell', {folder_files.name}, 'dictionary_folder_filename','','Mandatory',...
        {'Select the folder containing Pre/Post dico files (.json), ratio dico file (.mat) and/or model file (.mat)'}};
    
    user_parameter(:,5)   = {'   .Prefix','char', '', 'prefix', '', '',...
        {'Choose a prefix for your output maps'}};
    
    user_parameter(:,8)   = {'   .Smooth?','cell', {'Yes','No'}, 'filtered', '', '',...
        {'Select ''Yes'' to smooth the signals  (recommanded ''No'')'}};
    user_parameter(:,9)   = {'   .Remove last echoes?','numeric','', 'removed', '', '',...
        {''}};
    
    user_parameter(:,6)   = {'   .Parameters','check', ...
        {'Vf', 'VSI', 'R', 'SO2', 'DH2O', 'B0theta', 'khi', 'Hct', 'T2', 'B0off', 'Label'},...
        'Params', '', '',...
        {'Select the parameters considered in the model'
        }'};
    user_parameter(:,7)   = {'   .Clip','char', '0 inf' 'clipValues', '', '',...
        {'Clip produced maps to desired values. Indicate "min(param 1), max(param 1); min(param 2), max(param 2); ...", default is no clip'}};
    user_parameter(:,10)   = {'   .Confidence Maps','cell', {'Yes','No'}, 'confMaps', '', '',...
        {'Produce confidence maps, default = No'}};
    user_parameter(:,11)   = {'   .Method','cell', {'DBL', 'DBM'}, 'method', '', '',...
        { 'Select:'
        '	- ''DBL'' to use the dictionary-based learning reconstruction method [Fabien Boux et al.]'
        '	- ''DBM'' to use the dictionary-based matching reconstruction method [Dan Ma et al.]'
        }'};
    user_parameter(:,12)   = {'   .Make ScoreMap?','cell', {'Yes','No'}, 'mkScoreMap', '', '',...
        {'Select ''Yes'' add score map in outputs  (recommanded ''Yes'')'; 'WARNING: does not work with DBL'}'};
    user_parameter(:,13)  = {'   .Model settings (if the DBL method is used):','Text','','','','',...
        {'Recommanded:'
        'K = 50'
        'Lw = 0'
        'cstr on Sigma = ''d*'''
        'cstr on Gamma = ''  '''
        'Dictionary augmentation = 60'
        }'};
    user_parameter(:,14)  = {'       .Number of regions','numeric','','K','','',...
        {'Recommanded: K = 50'
        'If K is -1, an automatic tuning of the parameter is performed using BIC (time-consuming)'
        }'};
    user_parameter(:,15)  = {'       .Number of additional unsupervised parameter','numeric','','Lw','','',...
        {'Recommanded: Lw = 0'
        'If Lw is -1, an automatic tuning of the parameter is performed using BIC (time-consuming)'
        }'};
    user_parameter(:,16)  = {'       .Model constraint on Sigma','cell',{'i*','i','d*','d',' '},'cstrS','','',...
        {'''d'' = diagonal'
        '''i'' = isotropic'
        '''*'' = equal for all K regions'
        }'};
    user_parameter(:,17)  = {'       .Model constraint on Gamma','cell',{'i*','i','d*','d',' '},'cstrG','','',...
        {'''d'' = diagonal'
        '''i'' = isotropic'
        '''*'' = equal for all K regions'
        }'};
    user_parameter(:,18)   = {'       .Dictionary augmentation?','numeric','', 'augment', '', '',...
        {'Recommanded: SNR = 60'
         'This value represents the signal-to-noise ratio (SNR) between dictionary signals and noise added on these signals'
         'Inf value leads to no augmentation by noise addition (not recommanded)'
         }};
    
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
    
    if strcmp(opt.method,'DBL') && strcmp(opt.confMaps,'Yes')
        nb = i;
        
        for i = 1:numel(opt.Params)
            opt.Table_out(nb+i,:) = opt.Table_out(1,:);

            opt.Table_out.Filename(nb+i) =  categorical(cellstr([char(opt.Table_out.Patient(i)), '_', char(opt.Table_out.Tp(i)), '_', char(opt.Table_out.SequenceName(i)) '_confidence'])); 
            opt.Table_out.SequenceName(nb+i) = categorical(cellstr([char(opt.prefix), char(opt.Params{i}), '_confidence']));

            files_out.In2{i} = [char(opt.Table_out.Path(nb+i)), char(opt.Table_out.Filename(nb+i)) '.nii'];
        end
    end
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
if ~isempty(opt.clipValues)
    clipValues = str2num(opt.clipValues);
    assert(numel(opt.Params) == size(clipValues, 1), 'The number of clipping couples (min, max) provided does not match the number of parameters')
    flagClip = 1;
else
    flagClip = 0;
end

% Read json files and create ratio dictionary
opt.dictionary_folder_filename = fullfile(opt.folder, opt.dictionary_folder_filename);

d   = dir(opt.dictionary_folder_filename);
opt.dictionary_pre_filename     = d(contains({d.name}, 'PRE_'));
if ~isempty(opt.dictionary_pre_filename)
    opt.dictionary_pre_filename     = opt.dictionary_pre_filename.name;
end
opt.dictionary_post_filename    = d(contains({d.name}, 'POST_'));
if ~isempty(opt.dictionary_post_filename)
    opt.dictionary_post_filename    = opt.dictionary_post_filename.name;
end


% Name the model using a random ID and check if this model already exist
if strcmp(opt.method, 'DBL')
    model_filename = [opt.dictionary_folder_filename filesep 'MODEL_' num2str(10^5+randi(10^6-10^5,1)) '.mat'];

    list_models = d(contains({d.name}, 'MODEL_'));
    flag_model_exist = exist(model_filename,'file');
    flag_valid = false(size(list_models));
    for i = 1:length(list_models)

        load([opt.dictionary_folder_filename filesep list_models(i).name],'Parameters')
        
        try
            flag_valid(i) = strcmp(opt.cstrG, Parameters.cstr.Gammat);
            flag_valid(i) = flag_valid(i) && strcmp(opt.cstrS, Parameters.cstr.Sigma);
            flag_valid(i) = flag_valid(i) && opt.K == Parameters.K;
            flag_valid(i) = flag_valid(i) && opt.Lw == Parameters.Lw;
            flag_valid(i) = flag_valid(i) && opt.augment == Parameters.data_augmentation;
        catch
            flag_valid(i) = 0;
        end
    end
    flag_model_exist = flag_model_exist || any(flag_valid);
    if flag_model_exist
        l = find(flag_valid == true);
        model_filename = [opt.dictionary_folder_filename filesep list_models(l(1)).name];
    end
end


% If no model has already be computed create or load dictionary
if (strcmp(opt.method, 'DBL') && ~flag_model_exist) || strcmp(opt.method, 'DBM')
    
    % Dictionary filename
    dico_filename = [opt.dictionary_folder_filename filesep 'DICO.mat'];
    
    % If, dico exists, load it, else, create it and save it
    if exist(dico_filename,'file')
        load(dico_filename)
        
    else
        Pre     = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_pre_filename]);
        Post    = loadjson([opt.dictionary_folder_filename filesep opt.dictionary_post_filename]);

        Dico.MRSignals      = abs(Post.MRSignals) ./ abs(Pre.MRSignals); % Ratio post/pre signals 
        Dico.Tacq           = Pre.Sequence.Tacq;
        Dico.Parameters.Par = Pre.Parameters.Par; % Parameters used to simulate X signals
        Dico.Parameters.Labels = Pre.Parameters.Labels;
        clear Pre Post
        save(dico_filename,'Dico')
    end
end


% Load and if necessary, smooth observations
Xpost   = niftiread(files_in.In2{1});
Xpre    = niftiread(files_in.In1{1});

% Generate ratio signals from scans (and TODO: ROI if given)
Xobs            = Xpost ./ Xpre;
json_filename   = split(files_in.In2{1},'.');
json_filename{2} = '.json';
Obs             = ReadJson([json_filename{1} json_filename{2}]);

% Smooth
if strcmp(opt.filtered, 'Yes') == 1
    sigma = .5;
    
    FilteredImages = nan(size(Xobs));
    for z = 1:size(Xobs,3); for t = 1:size(Xobs,4)
        FilteredImages(:,:,z,t) = imgaussfilt(Xobs(:,:,z,t), sigma, 'FilterSize',3);
    end; end
    Xobs = FilteredImages;
    clear FilteredImages
end


% If necessary, remove echoes
if opt.removed > 0
    clear tmp
    for x = 1:size(Xobs,1); for y = 1:size(Xobs,2); for z = 1:size(Xobs,3)
        tmp(x,y,z,:) = Xobs(x,y,z,1:end-opt.removed);
    end; end; end
    Xobs    = tmp;
    Obs.EchoTime.value = Obs.EchoTime.value(1:end-opt.removed);
end

% Permute
Xobs(Xobs == inf) = nan;
Xobs	= permute(Xobs, [1 2 4 3]);


% Reformat dico (not needed if MODEL is already computed)
if strcmp(opt.method, 'DBM') || ~flag_model_exist
    tmp = nan(size(Dico.MRSignals,1), length(Obs.EchoTime.value'));
    if size(Xobs,length(size(Xobs))) ~= size(Dico.MRSignals,2)
        warning('Sizes of scans and dictionary MR signals are differents: dictionary MR signals reshaped')
        for i = 1:size(Dico.MRSignals,1)
            tmp(i,:) = interp1(Dico.Tacq(1:size(Dico.MRSignals,2)), Dico.MRSignals(i,:), Obs.EchoTime.value'*1e-3);
        end
    end
    Dico.MRSignals = tmp;
    Dico.Tacq   = Obs.EchoTime.value'*1e-3;
    %remove row containning nan values
    nn = ~any(isnan(Dico.MRSignals),2);
    Dico.MRSignals  = Dico.MRSignals(nn,:);
    Dico.Parameters.Par = Dico.Parameters.Par(nn,:);
    TmpDico{1}      = Dico;
end


% Compute DBM/DBL method
switch opt.method
    
    case 'DBM'
        % TODO: find something nicer than this permute trick
        Estimation  = AnalyzeMRImages(Xobs,TmpDico,opt.method,[]);
        Map.Y       = permute(Estimation.GridSearch.Y, [1 2 4 3]);
        
    case 'DBL'
        
        backgroundROI = false(size(Xpre,1), size(Xpre,2));
        x = ceil(size(Xpre,1) .* [0.02 0.17]);
        y = ceil(size(Xpre,2) .* [0.02 0.17]);
        backgroundROI(x(1):x(2),y(1):y(2)) = true;
        
%         SNRmap = 1./ ( 1./computeSNRmap(squeeze(Xpost) ,backgroundROI) ...
%             + 1./ computeSNRmap(squeeze(Xpre),backgroundROI) );

%         SNRmap = 1./ ( 1./computeSNRmap(squeeze(permute(Xpost, [1 2 4 3])) ,backgroundROI) ...
%             + 1./ computeSNRmap(squeeze(permute(Xpre, [1 2 4 3])),backgroundROI) );

        SNRmap = [];
        % Compute the learning only one time per dictionary
        if exist(model_filename,'file')
            load(model_filename,'Parameters','labels');
            Estimation = AnalyzeMRImages(Xobs, [], opt.method, Parameters,[],[], SNRmap);
            
            TmpDico{1}.Parameters.Labels = labels;
        
        else
            count = 1;
            
            % Pick parameters of interest
            for i = 1:length(Dico.Parameters.Labels)
                tmp = split(Dico.Parameters.Labels{i},'.',2);
                if any(strcmp(tmp{end}, opt.Params))
                    params.Par(:,count)     = TmpDico{1}.Parameters.Par(:,i);
                    params.Labels{count}    = TmpDico{1}.Parameters.Labels{i};
                    count = count +1;
                end
            end
            TmpDico{1}.Parameters = params;
            
            
            if opt.augment ~= 0
                TmpDico{1}.MRSignals = AddNoise(TmpDico{1}.MRSignals, opt.augment);
            end
            
            % Train regression
            % Parameters of the regression
            clear Parameters
            if opt.K >= 0,  Parameters.K = opt.K; end
            if opt.Lw >= 0, Parameters.Lw = opt.Lw; end
            if strcmp(opt.cstrS,' '), opt.cstrS = 'd*'; end
            if strcmp(opt.cstrG,' '), opt.cstrG = ''; end
            Parameters.cstr.Sigma   = opt.cstrS;
            Parameters.cstr.Gammat  = opt.cstrG;
            Parameters.cstr.Gammaw  = '';
            Parameters.data_augmentation = opt.augment;
            Parameters.comb = 'Conc';
            
            [Estimation, Parameters] = AnalyzeMRImages(Xobs, TmpDico, opt.method, Parameters,[],[], SNRmap);
            
            labels      = TmpDico{1}.Parameters.Labels; 
            save(model_filename,'Parameters', 'labels')
        end
            
        Map.Y   = permute(Estimation.Regression.Y, [1 2 4 3]);
        Map.Std	= permute(Estimation.Regression.Cov, [1 2 4 3]).^.5;
end


% Extract maps (and modify unit if necessary)
count = 1;
for i = 1:length(TmpDico{1}.Parameters.Labels)
    tmp = split(TmpDico{1}.Parameters.Labels{i},'.',2);
    if any(strcmp(tmp{end}, opt.Params))
        
        Labels{count} = tmp{end};      
        idxParam = find(strcmp(tmp{end}, opt.Params));
        
        switch tmp{end}
        %If the ression method is performed, extract also the confidence maps
            case {'Vf', 'SO2'} %convert to percent 
                MapStruct{count}    = 100*Map.Y(:,:,:,i);
                if strcmp(opt.method, 'DBL') && strcmp(opt.confMaps,'Yes')
                    StdStruct{count}    = 100*Map.Std(:,:,:,i);
                end
                if flagClip
                    Min = clipValues(idxParam, 1);
                    Max = clipValues(idxParam, 2);
                    MapStruct{count}(MapStruct{count} < Min) = Min;
                    MapStruct{count}(MapStruct{count} > Max) = Max;
                end
            case {'VSI', 'R'} % convert m to µm
                MapStruct{count}    = 1e6*Map.Y(:,:,:,i);
                if strcmp(opt.method, 'DBL') && strcmp(opt.confMaps,'Yes')
                    StdStruct{count}	= 1e6*Map.Std(:,:,:,i);
                end
                if flagClip
                    Min = clipValues(idxParam, 1);
                    Max = clipValues(idxParam, 2);
                    MapStruct{count}(MapStruct{count} < Min) = Min;
                    MapStruct{count}(MapStruct{count} > Max) = Max;
                end
            case {'DH2O'} % convert 
                MapStruct{count}    = 1e12*Map.Y(:,:,:,i);
                if strcmp(opt.method, 'DBL') && strcmp(opt.confMaps,'Yes')
                    StdStruct{count}	= 1e12*Map.Std(:,:,:,i);
                end
                if flagClip
                    Min = clipValues(idxParam, 1);
                    Max = clipValues(idxParam, 2);
                    MapStruct{count}(MapStruct{count} < Min) = Min;
                    MapStruct{count}(MapStruct{count} > Max) = Max;
                end
            case {'T2'} % convert 
                MapStruct{count}    = 1e3*Map.Y(:,:,:,i);
                if strcmp(opt.method, 'DBL') && strcmp(opt.confMaps,'Yes')
                    StdStruct{count}	= 1e3*Map.Std(:,:,:,i);
                end
                if flagClip
                    Min = clipValues(idxParam, 1);
                    Max = clipValues(idxParam, 2);
                    MapStruct{count}(MapStruct{count} < Min) = Min;
                    MapStruct{count}(MapStruct{count} > Max) = Max;
                end
            otherwise
                MapStruct{count} 	= Map.Y(:,:,:,i);
                if strcmp(opt.method, 'DBL') && strcmp(opt.confMaps,'Yes')
                    StdStruct{count}    = Map.Std(:,:,:,i);
                end
                if flagClip
                    Min = clipValues(idxParam, 1);
                    Max = clipValues(idxParam, 2);
                    MapStruct{count}(MapStruct{count} < Min) = Min;
                    MapStruct{count}(MapStruct{count} > Max) = Max;
                end
        end
        
        count   = count +1;
    end
end

if strcmp(opt.mkScoreMap, 'Yes')
    MapStruct{end+1} = Estimation.scoreMap;
    Labels{end+1} = 'ScoreMap';
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
                
        if contains(files_out.In1{i}, [Labels{j}, '.nii'])
            [path, name, ~] = fileparts(files_out.In1{i});
            WriteJson(J, [path, '/', name, '.json'])
            
            info.Filename = files_out.In1{i};
            info.ImageSize = size(MapStruct{j});
            info.PixelDimensions = info.PixelDimensions(1:length(size(MapStruct{j})));
            info.Datatype = class(MapStruct{j});
            niftiwrite(MapStruct{j}, files_out.In1{i}, info);


            if strcmp(opt.method, 'DBL')
                [path, name, ~] = fileparts(files_out.In2{i});
                WriteJson(J, [path, '/', name, '.json'])
                
                info.Filename = files_out.In2{i};
                info.ImageSize = size(StdStruct{j});
                info.PixelDimensions = info.PixelDimensions(1:length(size(StdStruct{j})));
                info.Datatype = class(StdStruct{j});
                niftiwrite(StdStruct{j}, files_out.In2{i}, info);
            end
        end
    end
end