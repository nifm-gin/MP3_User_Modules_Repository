function [files_in,files_out,opt] = Module_Export_BIDS_like(files_in,files_out,opt)

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
    module_option(:,3)   = {'RefInput',1};
    module_option(:,4)   = {'InputToReshape',1};
    module_option(:,5)   = {'Table_in', table()};
    module_option(:,6)   = {'Table_out', table()};
    module_option(:,7)   = {'Output_Folder','Export_BIDS'};
    module_option(:,8)   = {'AutomaticJobsCreation', 'No'};
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
         {'This module copy every scan selected in a BIDS''s like format :'
         'Output_Folder'
         '    .csv containing the information of the database'
         '    --> sub_(patient''name)'
         '        --> Anat'
         '            --> sub_(patient''name)_(scan''s name).nii.gz'
         }
        };
    user_parameter(:,2)   = {'   .Scan','XScan','','', {'SequenceName'},'Mandatory',...
         'Please select the scans that will be exported'};
    user_parameter(:,3)   = {'   .Output folder','char','','Output_Folder','', '','the name of the folder where will be copied your data.'};
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional', 'Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)', user_parameter(7,:)', 'VariableNames', VariableNames);
%%
    
    % So for no input file is selected and therefore no output
    % The output file will be generated automatically when the input file
    % will be selected by the user
    files_in = {''};
    files_out = {''};
    return
  
end
%%%%%%%%

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Export_Values_VoxelByVoxel:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%create Output_Folder

if ~exist(strrep(opt.folder_out, '/Tmp', [filesep opt.Output_Folder]), 'dir')
    mkdir(strrep(opt.folder_out, '/Tmp', [filesep opt.Output_Folder]))
end
% save as cvs
cvs_table = table;
cvs_table.Group =  opt.Table_in.Group;
cvs_table.Patient =  opt.Table_in.Patient;
cvs_table = unique(cvs_table);
cvs_table.Patient = categorical(strcat('sub-', cellstr(cvs_table.Patient)));


f_out = strrep(opt.folder_out, '/Tmp', [filesep opt.Output_Folder, filesep, 'metadata.csv']);

writetable(cvs_table, f_out) 



patient_list = unique(opt.Table_in.Patient);
for i=1:numel(patient_list)
    anat_path = strrep(opt.folder_out, '/Tmp', [filesep opt.Output_Folder filesep 'sub-' char(patient_list(i)) filesep 'anat']);
    tmp_table = opt.Table_in(opt.Table_in.Patient == patient_list(i),:);
    for j = 1 : size(tmp_table,1)
        gzip(strcat(char(tmp_table.Path(j)), char(tmp_table.Filename(j)), '.nii'), anat_path)
        movefile(strcat(anat_path, filesep, char(tmp_table.Filename(j)), '.nii.gz'), strcat(anat_path, filesep, 'sub-', char(patient_list(i)), '_', char(tmp_table.SequenceName(j)), '.nii.gz'))
    end
end

