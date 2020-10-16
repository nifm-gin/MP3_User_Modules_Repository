function [files_in,files_out,opt] = Module_SPM_old_normaliseWrite(files_in,files_out,opt)

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
    module_option(:,3)   = {'Preserve','Preserve Concentrations'};
    module_option(:,4)   = {'Bounding_box','-78 -112  -70; 78   76   85'};
    module_option(:,5)   = {'Voxel_sizes','2 2 2'};
    module_option(:,6)   = {'Interpolation','Trilinear'};
    module_option(:,7)   = {'Wrapping','No wrap'};
    module_option(:,8)   = {'output_filename_prefix','w'};
    module_option(:,9)   = {'RefInput',1};
    module_option(:,10)   = {'InputToReshape',1};
    module_option(:,11)   = {'Table_in', table()};
    module_option(:,12)   = {'Table_out', table()};
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
    user_parameter(:,1)   = {'Description','Text','','','', '','Write out warped images : refered as old normalisation write in spm12'}  ;
    user_parameter(:,2)   = {'Images to transform','XScanOrXROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Transformation information','1mat','','',{'SequenceName'}, 'Mandatory','the transformation is store in a _seg_inv_sn.mat or _seg_sn.mat file obtained from the module SPM_old_segment'};
    user_parameter(:,4)   = {'Writing options','','','','','',''};
    user_parameter(:,5)   = {'   .Preserve','cell',{'Preserve Concentrations', 'Preserve Amount'},'Preserve','','',...
        {'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.'
         ''
         'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity.'
        }'};
    user_parameter(:,6)   = {'   .Bounding box','numeric','','Bounding_box','','',...
        {'bounding box (2x3 matrix - in mm)'
         'Non-finite values mean use template bounding box (user can enter NaN NaN NaN; NaN NaN NaN to use the template bounding box)' 
         'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'
        }'};
    user_parameter(:,7)   = {'   .Voxel sizes','numeric','','Voxel_sizes','', '',...
        {'The voxel sizes (x, y & z, in mm) of the written normalised images.'
        'Non-finite values mean use template vox (user can enter NaN NaN NaN to use the template voxel sizes)'}
        };
    user_parameter(:,8)   = {'   .Interpolation','cell',{'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'},'Interpolation','','',...
        'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation is slower because it uses more neighbours.'};
    user_parameter(:,9)   = {'   .Warping','cell',{'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'},'Wrapping','','',...
        'These are typically: No wrapping - for PET or images that have already been spatially transformed.  Wrap in  Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'};
   
    user_parameter(:,10)   = {'   .Filename Prefix','char','w','output_filename_prefix','','',...
        'Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ''w''.'};
  
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
% select only the images to transform (the last row corresponds to the
% _sn.mat file
table_in_scans = opt.Table_in(1:end-1,:);


if isempty(files_out)
    opt.Table_out = table();
    for i = 1 : size(table_in_scans,1)
        New_instance_table_out = table_in_scans(i,:);
        New_instance_table_out.IsRaw = categorical(0);
        New_instance_table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
        New_instance_table_out.SequenceName = categorical(cellstr([opt.output_filename_prefix, char(New_instance_table_out.SequenceName)]));
        New_instance_table_out.Filename = categorical(cellstr([char(New_instance_table_out.Patient), '_', char(New_instance_table_out.Tp), '_', char(New_instance_table_out.SequenceName)]));
        f_out = [char(New_instance_table_out.Path), char(New_instance_table_out.Patient), '_', char(New_instance_table_out.Tp), '_', char(New_instance_table_out.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out ; New_instance_table_out];
    end
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_SPM_old_normaliseWrite:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end


%% If the test flag is true, stop here !

if opt.flag_test == 1
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The core of the brick starts here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Status, Message, Wrong_File] = Check_files(files_in);
if ~Status
    error('Problem with the input file : %s \n%s', Wrong_File, Message)
end


%  interpolation method (0-7)
switch opt.Interpolation
    case 'Nearest neighbour'
        flags.interp = 0;
    case 'Trilinear'
        flags.interp = 1;
    case '2nd Degree B-Spline'
        flags.interp = 2;
    case '3rd Degree B-Spline'
        flags.interp = 3;
    case '4th Degree B-Spline'
        flags.interp = 4;
    case '5th Degree B-Spline'
        flags.interp = 5;
    case '6th Degree B-Spline'
        flags.interp = 6;
    case '7th Degree B-Spline'
        flags.interp = 7;
end
% wrap edges (e.g., [1 1 0] for 2D MRI sequences)
switch opt.Wrapping
    case 'No wrap'
        flags.wrap = [0 0 0];
    case 'Warp X'
        flags.wrap = [1 0 0];
    case 'Warp Y'
        flags.wrap = [0 1 0];
    case 'Warp X&Y'
        flags.wrap = [1 1 0];
    case 'Warp Z'
        flags.wrap = [0 0 1];
    case 'Warp X&Z'
        flags.wrap = [1 0 1];
    case 'Warp Y&Z'
        flags.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        flags.wrap = [1 1 1];
end
% voxel sizes (3 element vector - in mm)
flags.vox = str2num(opt.Voxel_sizes); %#ok<ST2NM>
% bounding box (2x3 matrix - in mm)
flags.bb  = str2num(opt.Bounding_box); %#ok<ST2NM>
% either 0 or 1.  A value of 1 will "modulate"
if strcmp(opt.Preserve, 'Preserve Amount')
    flags.preserve = 1;
else
    flags.preserve = 0;
end
% Prefix for normalised images. Defaults to 'w'.
flags.prefix = opt.output_filename_prefix;


for i = 1:size(files_in.In1,2)
    % % First duplicate the reference image to the tmp folder and update the
    % files_in.In1 variable
    new_pathname = strcat(char(opt.Table_out.Path(1)),   char(opt.Table_in.Filename(i)), '.nii');
    if ~strcmp(files_in.In1{i},   new_pathname)
        copyfile(files_in.In1{i},   new_pathname)
        if exist(strrep(files_in.In1{i},'.nii','.json'), 'file')
            copyfile(strrep(files_in.In1{i},'.nii','.json'),  strrep(new_pathname,'.nii','.json'))
        end
    end
    
    files_in.In1{i} = new_pathname;
    
    VO = spm_write_sn(files_in.In1{i},files_in.In2{1},flags);
    VO.fname = files_out.In1{i};
    spm_write_vol(VO, VO.dat);
    %% Json Processing
    if exist(strrep(files_in.In1{i},'.nii','.json'), 'file')
        [path, name, ~] = fileparts(files_in.In1{i});
        jsonfile = [path, '/', name, '.json'];
        J = ReadJson(jsonfile);
        
        J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
        
        [path, name, ~] = fileparts(files_out.In1{i});
        jsonfile = [path, '/', name, '.json'];
        WriteJson(J, jsonfile)
    end
    
end
