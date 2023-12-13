function [files_in,files_out,opt] = Module_MergeScoreMaps(files_in,files_out,opt)

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
    module_option(:,4)   = {'Constant',1};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    module_option(:,9)   = {'output_filename_ext','_'};

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
    user_parameter(:,1)   = {'Description','Text','','','', '',...
        {'This module takes as inputs different scoremaps from several matches' 
        'It compute the max of all of these scoremaps and return an index map to know in which scoremap is the max for each voxel.'
        'The index map could then be used with the module "MergeMapsFromIndexes'
        'Be careful that the names of the score maps are sorted correctly as the parameters maps.'}' } ;
    user_parameter(:,2)   = {'Select the first scoremap','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the other scans','XScan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Smooth','output_filename_ext','','',...
        {'Specify the string to be added to the output filename.'
        'Default filename is ''ScoreMap_Indexes_Merged''.'}'};
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

if isempty(files_out)
    opt.Table_out = opt.Table_in(1,:);
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr(['ScoreMap_Indexes_Merged', opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Arithmetic:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
n_maps = size(files_in.In2,2) + 1 ;
for i=2:n_maps
    input(i).nifti_header = spm_vol(files_in.In2{i-1});
end

ref_scan = 1;

for i=1:n_maps
    Scoremap{i} = read_volume(input(i).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
%     Scoremap{i} = double(Scoremap{i});
end

%% max of score maps and index map
maxMatrix = Scoremap{1};
maxIndex = NaN(size(maxMatrix));
for i = 1:n_maps
    currentMatrix = Scoremap{i};
    % Find the maximum values among the matrices
    maxMatrix = max(maxMatrix, currentMatrix, 'omitnan');
    % Update the index matrix where the current matrix has a higher value
    maxIndex(currentMatrix == maxMatrix) = i;
end

% maxMatrix now contains the maximum values among the 6 matrices
% maxIndex contains the index of the matrix where the maximum value is located


% transform the maxIndex matrix in order to match to the nii header of the
% first input (rotation/translation)

maxIndex_reoriented = write_volume(maxIndex, input(1).nifti_header, 'Axial');

info = niftiinfo(files_in.In1{1});

nifti_header_output = info;
nifti_header_output.Filename = files_out.In1{1};
nifti_header_output.Filemoddate = char(datetime('now'));
[maxIndex_reoriented, FinalMat] = CropNifti(maxIndex_reoriented, nifti_header_output.Transform.T');
nifti_header_output.Datatype = class(maxIndex_reoriented);
nifti_header_output.Transform = affine3d(FinalMat');
nifti_header_output.ImageSize = size(maxIndex_reoriented); 
nifti_header_output.PixelDimensions = info.PixelDimensions(1:length(nifti_header_output.ImageSize));
nifti_header_output.MultiplicativeScaling = 1;

% % save the new .nii file
niftiwrite(maxIndex_reoriented, files_out.In1{1}, nifti_header_output);

%% Json processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
if isfile(jsonfile)
    J = ReadJson(jsonfile);
    
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    
    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end
