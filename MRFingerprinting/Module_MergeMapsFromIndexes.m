function [files_in,files_out,opt] = Module_MergeMapsFromIndexes(files_in,files_out,opt)

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
    module_option(:,9)   = {'output_filename_ext','Param'};

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
    user_parameter(:,1)   = {'Description','Text','','','', '', ...
        { 'This module merged different parameters maps into 1 based on an MaxIndex map' 
        'issued from ths scoremaps of all previous matches by the module MergeScoreMaps.'
        'Use this module for one parameters at time, and do not forget to specify the parameter in the extension name.'}' } ;
    user_parameter(:,2)   = {'Select indexed scoremap','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select parameters map','XScan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Smooth','output_filename_ext','','',...
        {'Specify the string to be added to the output filename.'
        'Default filename extension is ''_Merged''.'}'};
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
        opt.Table_out.SequenceName = categorical(cellstr(['Merged_match_', opt.output_filename_ext]));
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
niIndex.nifti_header = spm_vol(files_in.In1{1});

maxIndex = read_volume(niIndex.nifti_header,niIndex.nifti_header, 0, 'Axial');

n_maps = size(files_in.In2,2);

for i=1:n_maps
    input(i).nifti_header = spm_vol(files_in.In2{i});
end

ref_scan = 1;

for i=1:n_maps
    Maps{i} = read_volume(input(i).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
end

%% merged parameters maps according to index matrix

% Initialize FinalMap with zeros of the same size as the matrices in Maps
FinalMap = Maps{1};

% Iterate through each voxel and assign the corresponding value from the Maps cell array
for i = 1:size(FinalMap, 1)
    for j = 1:size(FinalMap, 2)
        for k = 1:size(FinalMap, 3)
            if isnan(maxIndex(i, j, k)) %If value at this iteration is nan, don't consider it
                continue
            end
            index = maxIndex(i, j, k);
            FinalMap(i, j, k) = Maps{index}(i, j, k);
        end
    end
end

FinalMap_reoriented = write_volume(FinalMap, input(1).nifti_header, 'Axial');

info = niftiinfo(files_in.In1{1});

nifti_header_output = info;
nifti_header_output.Filename = files_out.In1{1};
nifti_header_output.Filemoddate = char(datetime('now'));
[FinalMap_reoriented, FinalMat] = CropNifti(FinalMap_reoriented, nifti_header_output.Transform.T');
nifti_header_output.Datatype = class(FinalMap_reoriented);
nifti_header_output.Transform = affine3d(FinalMat');
nifti_header_output.ImageSize = size(FinalMap_reoriented); 
nifti_header_output.PixelDimensions = info.PixelDimensions(1:length(nifti_header_output.ImageSize));
nifti_header_output.MultiplicativeScaling = 1;

% % save the new .nii file
niftiwrite(FinalMap_reoriented, files_out.In1{1}, nifti_header_output);

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
