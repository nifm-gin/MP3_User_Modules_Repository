function [files_in,files_out,opt] = Module_cineASL(files_in,files_out,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the module's parameters with default values
if isempty(opt)
    %
    %     %%   % define every option needed to run this module
    %     % --> module_option(1,:) = field names
    %     % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'OutputSequenceName','Extension'};
    module_option(:,4)   = {'T1blood',1.7};
    module_option(:,5)   = {'T1tissue',1.5};
    module_option(:,6)   = {'Lambda',0.95};
    module_option(:,7)   = {'Alpha',8};
    module_option(:,8)   = {'Beta',0.5};
    module_option(:,9)   = {'Start_Index_TAG',8};
    module_option(:,10)   = {'Start_Index_CTR',8};
    module_option(:,11)   = {'RefInput',1};
    module_option(:,12)  = {'InputToReshape',1};
    module_option(:,13)  = {'Table_in', table()};
    module_option(:,14)  = {'Table_out', table()};
    module_option(:,15)  = {'output_filename_ext',''};
    module_option(:,16)  = {'NameOutFiles',''};
    
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','', '',{'Description of the module: '
        'The objective of this module is to create myocardial blood flow (MBF) maps from CineASL sequences.'
        'Details on the method can be found in this publication :'
        'Troalen T, Capron T, Bernard M, Kober F.'
        'In vivo characterization of rodent cyclic myocardial perfusion variation at rest and during adenosine-induced stress using cine-ASL cardiovascular magnetic resonance.'
        'J Cardiovasc Magn Reson. 2014 Feb 18;16(1):18.'
        'doi: 10.1186/1532-429X-16-18.'
        'PMID: 24548535;'
        'PMCID: PMC3937054'
        ''
        'This module was coded by B. Lemasson (GIN; Grenoble; France) using the code and assistance of F. Kober (CRMBM Marseille; France).'}
        }  ;
    user_parameter(:,2)   = {'Select the CineAsl-Control scan','1Scan','','',{'SequenceName'}, 'Mandatory','Please select a 5d scan correspondonding the the CineASL-control'};
    user_parameter(:,3)   = {'Select the CineAsl-Tag scan','1Scan','','',{'SequenceName'}, 'Mandatory','Please select a 5d scan correspondonding the the CineASL-Tag'};

    user_parameter(:,4)   = {'   .Filename exention','char','MBF','output_filename_ext','','',...
        {'Specify the filename extention to add the the following names :'
        'Output filenames are : '
        '     - MBF'
        '     - CineASL_AvCTL'
        '     - CineASLAvTag'
        '     - CineAslDeltaMA'}'};
    user_parameter(:,5)   = {'   .T1 of the blood','numeric', '','T1blood','','',{'Please entre the T1 value of the blood; defauld value is 1.7 sec'}'};
    user_parameter(:,6)   = {'   .T1 of the tissue','numeric', '','T1tissue','','',{'Please entre the T1 value of the tissue; defauld value is 1.5 sec'}'};
    user_parameter(:,7)   = {'   .Lambda','numeric', '','Lambda','','',{'Please entre the Labda; defauld value is 0.95'}'};
    user_parameter(:,8)   = {'   .Alpha','numeric', '','Alpha','','',{'Please entre the alpha value; defauld value is 8'}'};
    user_parameter(:,9)   = {'   .Beta','numeric', '','Beta','','',{'Please entre the beta value; defauld value is 0.5'}'};
    user_parameter(:,10)   = {'   .Average start index (TAG)','numeric', '','Start_Index_TAG','','',{'Please enter the repetition from which the data will be analyzed for the Tagging'
        '(all images acquired before will be rejected).'}'};
    user_parameter(:,11)   = {'   .Average start index (Ctr)','numeric', '','Start_Index_CTR','','',{'Please enter the repetition from which the data will be analyzed for the Control'
        '(all images acquired before will be rejected).'}'};

    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
    %%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
    
end
%%%%%%%%
opt.NameOutFiles = {'MBF', 'CineAsl_DeltaMap', 'CineAsl_AvTag', 'CineAsl_AvCTL'};

if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in(opt.RefInput,:);
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
                table_out_tmp.SequenceName = categorical(cellstr([opt.NameOutFiles{i}, opt.output_filename_ext]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
end

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('ModuleCineASL:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%% load input Nii file

input(1).nifti_header = spm_vol(files_in.In1{1});
input(2).nifti_header =  spm_vol(files_in.In2{1});
cineASL_CTR = read_volume(input(1).nifti_header, input(opt.RefInput).nifti_header, 0);
cineASL_TAG = read_volume(input(2).nifti_header, input(opt.RefInput).nifti_header, 0);


input(1).json = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
%input(2).json = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));
format long
% Get the repetition time in sec
TR = input(1).json.RepetitionTime.value * 1e-3; 

% Get the flip angle value and convert it to radia
alpharad = deg2rad(input(1).json.FlipAngle.value);

% compute the T1star (which correspond to the T1 value under FLASH sequence)
T1star = opt.T1tissue * TR / (TR - opt.T1tissue * log(cos(alpharad))) ;

% calculate the steady state magnetization under FLASH sequence
mss = (1 - (exp(-TR/opt.T1tissue))) ./ (1 - (cos(alpharad) .* exp(-TR/opt.T1tissue)));

% Determine indexes to retain and to discard from marked selection
cineASL_CTR = cineASL_CTR(:,:,:,:,opt.Start_Index_CTR:end);
cineASL_TAG = cineASL_TAG(:,:,:,:,opt.Start_Index_TAG:end);

% Average over cardiac cycles
cineASL_CTR_av  = mean(cineASL_CTR,  5);
cineASL_TAG_av  = mean(cineASL_TAG,  5);

%compute the absolute and the relative difference images
absdelta = cineASL_CTR_av - cineASL_TAG_av;
delta = absdelta ./ cineASL_CTR_av .* (cineASL_TAG_av ~= 0);


%calculate the MBF according to Capron et al. MRM 2013
MBF = opt.Lambda * mss ./ T1star * delta ./ (2. * opt.Beta - delta) * 60;       


% Map to save :  MBF, delta, series_AvTag, series_AvCtrl
mapsVar = {MBF, delta, cineASL_TAG_av, cineASL_CTR_av};
info  =  niftiinfo(files_in.In1{1});
for i=1:length(mapsVar)
    Images_reoriented = write_volume(mapsVar{i}, input(opt.RefInput).nifti_header, 'axial');

    info_to_save =  info;
    info_to_save.Filename = files_out.In1{i};
    info_to_save.Filemoddate = char(datetime('now'));
    info_to_save.Datatype = class(Images_reoriented);
    info_to_save.PixelDimensions = info.PixelDimensions(1:length(size(Images_reoriented)));
    info_to_save.ImageSize = size(Images_reoriented);



    niftiwrite(Images_reoriented, files_out.In1{i}, info_to_save)

    %% Json Processing
    new_json = KeepModuleHistory(input(1).json , struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
    [path, name, ~] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(new_json, jsonfile)

end


%


