function [files_in,files_out,opt] = Module_permute_4th_and_3rd_dim(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename_ext','_permuted'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'Slice_thickness_option','TR'};
    module_option(:,6)   = {'Slice_thicness_manually_added', ' '};
    module_option(:,7)   = {'RefInput',1};
    module_option(:,8)   = {'InputToReshape',1};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)   = {'Table_out', table()};
    
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
    user_parameter(:,1)   = {'Description','Text','','','', '', {'This module permute the 4th dimension with the 3rd dimension of a 4D scan'
                                                                   'the user must add the thickness of the slice that the new 3rd dimension will have (formerly the 4th)'
                                                                   'the new slice thickness could be automatically determined to correspond to the : '
                                                                   '     - Repetition Time (TR)'
                                                                   '     - Echo Time (TE)'
                                                                   '     - Inversion Time'
                                                                   '     - or manually added'
                                                                   ''
                                                                   'WARNING !!!!'
                                                                   'MP3 can not handle scans with variable thickness, therefore if you use the TR values as the new slice thicknesses (or TE or inversion time)'
                                                                   'this module will compute and use the mean TR value in between 2 images'}};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Permuted','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Select the new slice thickness','cell', {'TR', 'TE', 'Inversion Time', 'Manually added'},'Slice_thickness_option',...
        '', '','If Manual value is selected, please enter its value below'};
    user_parameter(:,6)   = {'        .Slice thickness (if manually added is selected)','numeric', '', 'Slice_thickness_manually_added', '', '', ''};

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
    error('permute_4th_and_3rd_dim:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get data from a specific axe      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info_spm = spm_vol(files_in.In1{1});
N = read_volume(info_spm, info_spm, 0);


info = niftiinfo(files_in.In1{1});

[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
if isfile(jsonfile)
    J = ReadJson(jsonfile);
end

% determine the new slice thickness
switch opt.Slice_thickness_option
    case 'TR'
         if length(J.RepetitionTime.value) > 1
            new_slice_thickness = (max(J.RepetitionTime.value)-min(J.RepetitionTime.value))/(length(J.RepetitionTime.value)-1);
         else
            new_slice_thickness = double(J.RepetitionTime.value);
         end
    case 'TE'
        if length(J.EchoTime.value) > 1
            new_slice_thickness = (max(J.EchoTime.value)-min(J.EchoTime.value))/(length(J.EchoTime.value)-1);
        else
            new_slice_thickness = J.EchoTime.value;
        end
       
    case 'Inversion Time'
        % inversion time
        if ~isfield(J, 'FairTIRArr')
            new_slice_thickness = max(J.FairTIRArrPVM.value)/length(J.FairTIRArrPVM.value);
        else
            new_slice_thickness = max(J.FairTIRArr.value)/length(J.FairTIRArr.value);
        end
        
    case 'Manually added'
        if strcmp(opt.Slice_thicness_manually_added, ' ')
            error('No slice thickness value has been filled : %s \n%s', Wrong_File, Message)
        else
            new_slice_thickness = str2double(opt.Slice_thicness_manually_added);
        end
end


info_permuted = info;
% if only one slice was acquired, we have to determine the where the slice is located (which dimension) 
dimension_of_the_slices = find(size(N) ==1);
if isempty(dimension_of_the_slices)
    dimension_of_the_slices = 0;
end
switch dimension_of_the_slices
    case 1 % the slice is located in the first dimension
        N_permuted = permute(N, [2 3 4 1]);
        N_permuted_reoriented = flip(N_permuted,1);
        N_permuted_reoriented = flip(N_permuted_reoriented,3);

        info_permuted.ImageSize = size(N_permuted_reoriented);
        info_permuted.PixelDimensions(4) = info_permuted.PixelDimensions(3);
        info_permuted.PixelDimensions(3) = new_slice_thickness;
        info_permuted.PixelDimensions = info_permuted.PixelDimensions(1:length(size(N_permuted_reoriented)));
    case 2 % the slice is located in the second dimension
        N_permuted_reoriented = permute(N, [1 3 4 2]);
        N_permuted_reoriented = flip(N_permuted_reoriented,3);

        info_permuted.ImageSize = size(N_permuted_reoriented);
        info_permuted.PixelDimensions(4) = info_permuted.PixelDimensions(3);
        info_permuted.PixelDimensions(3) = new_slice_thickness;
        info_permuted.PixelDimensions = info_permuted.PixelDimensions(1:length(size(N_permuted_reoriented)));        
      
    otherwise  % the slice is located in the third dimension
        N_permuted = permute(N, [1 2 4 3]);
        N_permuted_reoriented = write_volume(N_permuted, info_spm, 0);
        
        info_permuted.ImageSize = size(N_permuted_reoriented);
        info_permuted.PixelDimensions(4) = info_permuted.PixelDimensions(3);
        info_permuted.PixelDimensions(3) = new_slice_thickness;
        info_permuted.PixelDimensions = info_permuted.PixelDimensions(1:length(size(N_permuted_reoriented)));
        
end
      

% complet information saved
info_permuted.Filemoddate = char(datetime('now'));
info_permuted.Description = [info_permuted.Description, 'Modified by Module_permute_4th_and_3rd_dim'];
    


niftiwrite(N_permuted_reoriented, files_out.In1{1} , info_permuted)

%%% work to do : update the json in accordance to the permuation made 
%%% for instance EchoTime should be seted to 0 ??
if exist(jsonfile)
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end
