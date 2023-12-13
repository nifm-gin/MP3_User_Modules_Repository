function [files_in,files_out,opt] = Module_ProcessBadOSSI(files_in,files_out,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Initialize the module's parameters with default values 
if isempty(opt)

%     %   % define every option needed to run this module
%     % --> module_option(1,:) = field names
%     % --> module_option(2,:) = defaults values    
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'output_filename_ext','_Corrected'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'SliceNumber', 1};
    module_option(:,6)   = {'RefInput',1};
    module_option(:,7)   = {'InputToReshape',1};
    module_option(:,8)   = {'Table_in', table()};
    module_option(:,9)   = {'Table_out', table()};
    
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
% 
        %% list of everything displayed to the user associated to their 'type'
         % --> user_parameter(1,:) = user_parameter_list
         % --> user_parameter(2,:) = user_parameter_type
         % --> user_parameter(3,:) = parameter_default
         % --> user_parameter(4,:) = psom_parameter_list
         % --> user_parameter(5,:) = Scans_input_DOF : Degrees of Freedom for the user to choose the scan
         % --> user_parameter(6,:) = IsInputMandatoryOrOptional : If none, the input is set as Optional. 
         % --> user_parameter(7,:) = Help : text data which describe the parameter (it
         % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','', '', {'Removes unwanted null echoes from poorly written OSSI sequence'}};
    user_parameter(:,2)   = {'Select one scan as input','1ScanOr1ROI','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_S2E','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Slice Number', 'numeric', 1, 'SliceNumber','', '',{'Number of slices in the volume'}};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)','VariableNames', VariableNames);
%%
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
  
end
%%%%%%%%

if isempty(files_out)
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);   
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.output_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Resphape:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

% [path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
% if ~strcmp(ext_nii, '.nii')
%      error('First file need to be a .nii, not a %s. Path : %s, Name : %s, opt.files_in : %s', ext_nii, path_nii, ext_nii, files_in);  
% end
% 
% %% Building default output names
% if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
%     opt.folder_out = path_nii;    
% end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get data from a specific axe      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open the data in the correct orientation (as the display in MP3
info_spm = spm_vol(files_in.In1{1});
N = read_volume(info_spm, info_spm, 0);

% nii hearder
info = niftiinfo(files_in.In1{1});

% N = niftiread(files_in.In1{1});
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];


Informations = whos('N');
% axes = str2double(opt.Dimenssion_to_reduce);
% index = str2num(opt.Index_to_keep);

[X, Y, nImages] = size(N);
NewIm = cat(3, N(:,:,1:150), N(:,:,251:400), N(:,:,501:650), N(:,:,751:900));
% reorient data before saving it
NewIm_reoriented = write_volume(NewIm, info_spm, 0);

info2 = info;

NewIm = cast(NewIm_reoriented,info.Datatype); 
info2.ImageSize = size(NewIm_reoriented);
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Description = [info.Description, 'Modified by Slice2Echo module'];


niftiwrite(NewIm_reoriented, files_out.In1{1}, info2)

if exist(jsonfile)
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    if axes == 4 && ~contains(J.ProtocolName.value, 'ADC') && numel(J.EchoTime.value) > numel(index(1):index(end))
        J.EchoTime.value = J.EchoTime.value(index(1):index(end));
    end

    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end
