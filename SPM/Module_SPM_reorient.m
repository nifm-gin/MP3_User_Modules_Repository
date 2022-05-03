function [files_in,files_out,opt] = Module_SPM_reorient(files_in,files_out,opt)

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
    module_option(:,3)   = {'right',0};
    module_option(:,4)   = {'forward',0};
    module_option(:,5)   = {'up',0};
    module_option(:,6)   = {'pitch',0};
    module_option(:,7)   = {'roll',0};
    module_option(:,8)   = {'yaw',0};
    module_option(:,9)   = {'resize_x',1};
    module_option(:,10)   = {'resize_y',1};
    module_option(:,11)   = {'resize_z',1};
    module_option(:,12)   = {'output_filename_ext','_Reoriented'};
    module_option(:,13)   = {'RefInput',1};
    module_option(:,14)   = {'InputToReshape',1};

    
    module_option(:,15)   = {'Table_in', table()};
    module_option(:,16)   = {'Table_out', table()};
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
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {' This module allows to reorient an image.'
        ''
        'The image can be re-oriented by entering appropriate translations, rotations and zooms into the panel on the left.'
        'The transformations are considered to be relative to any existing transformations that may be stored.'
        ''
        'For more information, please refere spm functions'}'};
    user_parameter(:,2)   = {'Reference Image','1ScanOr1ROI','','', {'SequenceName'},'Mandatory',...
        'Select the scans to modify'};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Reoriented','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .right (mm)','numeric',0,'right','', '',''};
    user_parameter(:,6)   = {'   .forward (mm)','numeric',0,'forward','', '',''};
    user_parameter(:,7)   = {'   .up (mm)','numeric',0,'up','', '',''};
    user_parameter(:,8)   = {'   .pitch (rad)','numeric',0,'pitch','', '',''};
    user_parameter(:,9)   = {'   .roll (rad)','numeric',0,'roll','', '',''};
    user_parameter(:,10)  = {'   .yaw (rad)','numeric',0,'yaw','', '',''};
    user_parameter(:,11)  = {'   .resize (x)','numeric',1,'resize_x','', '',''};
    user_parameter(:,12)  = {'   .resize (y)','numeric',1,'resize_y','', '',''};
    user_parameter(:,13)  = {'   .resize (z)','numeric',1,'resize_z','', '',''};

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
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
end



%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_SPM_reorient:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

% First duplicate the source scan using the prefix string (user-defined)
% Otherwise the spm will overwrite the file!!
copyfile(files_in.In1{1},  files_out.In1{1})

userdata  = [opt.right opt.forward opt.up opt.pitch opt.roll opt.yaw opt.resize_x opt.resize_y opt.resize_z 0 0 0];

matrice_transformation = spm_matrix(userdata);
 if det(matrice_transformation)<=0
        spm('alert!','This will flip the images',mfilename,0,1);
        return
 end
input_matrice = spm_get_space(files_out.In1{1});
spm_get_space(files_out.In1{1},matrice_transformation*input_matrice);


% update the json if needed
if exist(strrep(files_in.In1{1},'.nii','.json'), 'file')
    copyfile(strrep(files_in.In1{1},'.nii','.json'),  strrep(files_out.In1{1},'.nii','.json'))
    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    
    WriteJson(J, jsonfile) 
end
