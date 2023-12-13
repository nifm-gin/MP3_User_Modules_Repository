function [files_in,files_out,opt] = Module_iDREAM(files_in,files_out,opt)

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
    module_option(:,3)   = {'B1_filename_ext','_B1'};
    module_option(:,4)   = {'OutputSequenceName','Extension'};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    module_option(:,9)   = {'B1phi_filename_ext','_B1phase'};
    
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
    user_parameter(:,1)   = {'Description','Text','','','', '', {'Generates B1 and phase maps from Philips iDREAM standard acq.'}};
    user_parameter(:,2)   = {'Select one scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '', ''};
    user_parameter(:,4)   = {'   .B1 output filename extension','char','_B1','B1_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .df output filename extension','char','_df','B1phi_filename_ext','', '',''};
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
        opt.Table_out.SequenceName = categorical(cellstr(opt.B1_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.B1_filename_ext]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    
    dfMapRow.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        dfMapRow.SequenceName = categorical(cellstr(opt.B1phi_filename_ext));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        char(opt.Table_out.SequenceName)
        dfMapRow.SequenceName = categorical(cellstr([char(opt.Table_in.SequenceName), opt.B1phi_filename_ext]));
    end
    dfMapRow.Filename = categorical(cellstr([char(opt.Table_out.Patient(1)), '_', char(opt.Table_out.Tp(1)), '_', char(dfMapRow.SequenceName)]));
    dfMapRow.IsRaw = categorical(cellstr('0'));

    dfMapRow.Group = opt.Table_out.Group(1);
    dfMapRow.Patient = opt.Table_out.Patient(1);
    dfMapRow.Tp = opt.Table_out.Tp(1);
    dfMapRow.Type = opt.Table_out.Type(1);
    opt.Table_out = [opt.Table_out; struct2table(dfMapRow)];

    files_out.In1{end+1} = [char(opt.Table_out.Path(end)), char(opt.Table_out.Filename(end)) '.nii'];    
end


%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Resphape:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

if nImages == 80 % Philips calculated images
    B1vol = squeeze(N(:,:,1,:));
    dfvol = squeeze(N(:,:,2,:));
elseif nImages == 240 % Dicom raw + calculated images
    count = 1;
    B1vol = zeros(240,240,40);
    dfvol = zeros(240,240,40);
    for i = 1:40
        for j = 1:6
%             count = count +1;
%             count = (j-1) * 40 + i;
            if count <= 160  
                count = count+1;
                continue
            elseif count <= 200
                B1vol(:,:, count - 160) = N(:,:,i,j);
            else
                dfvol(:,:, count - 200) = N(:,:,i,j);    
            end
            count = count + 1;
        end 
    end 
    B1vol = cast(B1vol, 'single');
    dfvol = cast(dfvol, 'single');
end

B1vol = write_volume(B1vol, info_spm, 0);
dfvol = write_volume(dfvol, info_spm, 0);
    
info2 = info;
info2.PixelDimensions = info2.PixelDimensions(1:3);
info2.ImageSize = size(B1vol);
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Description = [info.Description, 'Modified by iDREAM module'];
info2.Filename = files_out.In1{1};

niftiwrite(B1vol, files_out.In1{1}, info2)
niftiwrite(dfvol, files_out.In1{2}, info2)

if exist(jsonfile)
    %B1
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    if axes == 4 && ~contains(J.ProtocolName.value, 'ADC') && numel(J.EchoTime.value) > numel(index(1):index(end))
        J.EchoTime.value = J.EchoTime.value(index(1):index(end));
    end

    [path, name, ~] = fileparts(files_out.In1{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
    
    %df
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    if axes == 4 && ~contains(J.ProtocolName.value, 'ADC') && numel(J.EchoTime.value) > numel(index(1):index(end))
        J.EchoTime.value = J.EchoTime.value(index(1):index(end));
    end

    [path, name, ~] = fileparts(files_out.In1{2});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end
