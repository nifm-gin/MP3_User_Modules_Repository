function [files_in,files_out,opt] = Module_DataToCpx(files_in,files_out,opt)

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
    module_option(:,3)   = {'output_filename_ext','my_sequence_Complex'};
    module_option(:,4)   = {'OutputSequenceName','AllName'};
    module_option(:,9)   = {'Constant',1};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,9)   = {'Table_in', table()};
    module_option(:,10)   = {'Table_out', table()};
    module_option(:,8)   = {'Normalize', 'No'};
    module_option(:,7)   = {'echoes', '3'};
    
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
    user_parameter(:,1)   = {'Description','Text','','','', '', {'This module is used to extract a sub-stack from a nD scan'}};
    user_parameter(:,2)   = {'Module input scan (M)','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Phase input scan (P)','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,4)   = {'   .Output filename extension','char','_Reshape','output_filename_ext','', '',''};
    user_parameter(:,5)   = {'   .Norm phase to first echo?','cell',{'Yes','No'}, 'Normalize', '', '',...
        {'Yes = for each TR, normalize the echoes phases to the first echo phase to avoid phase accumulation issue in matching. No = classic complex making'}};
    user_parameter(:,6)   = {'   .Echoes number?','char','', 'echoes', '', '',...
        {'Just specified echoes number to compute phase normalization correctly'}};
    
    
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
    opt.Table_out = opt.Table_in(1,:);
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

input(1).nifti_header = spm_vol(files_in.In1{1}); %module
input(2).nifti_header = spm_vol(files_in.In2{1}); %phase

input1 = read_volume(input(1).nifti_header, input(1).nifti_header, 0, 'Axial'); %module
input2 = read_volume(input(2).nifti_header, input(1).nifti_header, 0, 'Axial'); %phase

% %% phase unwrap 1
% input2(isnan(input2)) = 0;
% [Nx, Ny, Nz, Nt] = size(input2);
% img = zeros(Nx, Ny, Nz, Nt);
% for i=1:Nz
%     for j=1:Nt
% %         img(:,:,i,j) = unwrap_phase(img_wrapped(:,:,i,j));
%         img(:,:,i,j) = abs(unwrapLaplacian(input2(:,:,i,j), [Nx Ny]));
%     end
% end
% input2 = img

% %% phase unwrap 2
% % loop over all voxels in img
% [Nx, Ny, Nz, Nt] = size(input2);
% img = zeros(Nx, Ny, Nz, Nt);
% 
% for i = 1:Nx
%     for j = 1:Ny
%         for k = 1:Nz
%             if isnan(input2(i, j, k, 2))
%                 continue
%             end
%         img(i, j, k, :) = unwrap(squeeze(input2(i, j, k, :)));
%         end
%     end
% end
% input2 = img
%%
complex_data = input1 .* exp(1i * input2);
matrixSize = size(complex_data);
echoes = str2num(opt.echoes);

%% Dividing the echoes phase in a pulse by the first echo, trying to avoid phase accumulation issue in matching
% if strcmp(opt.Normalize, 'Yes') == 1
%     % Iterate over each signal in the matrix
%     for i = 1:matrixSize(1)
%         for j = 1:matrixSize(2)
%             for k = 1:matrixSize(3)
%                 signal = squeeze(complex_data(i, j, k, :)); % Extract the signal
%                 % Iterate over the signal in groups of n_echoes
%                 for n = 1:echoes:matrixSize(4)-(echoes-1)
%                     % Get the first point of the current group
%                     firstPoint = signal(n);
%                     % Compute the phase of the first point
%                     firstPhase = angle(firstPoint);
%                     % Divide the phase of each point in the group by the phase of the first point
%                     for m = n:(n + echoes-1)
%                         signal(m) = abs(signal(m)) * exp(1i * (angle(signal(m)) - firstPhase));
%                     end
%                 end
%                 % Store the modified signal back into the matrix
%                 complex_data(i, j, k, :) = signal;
%             end
%         end
%     end
% end

% Dividing the echoes phase on the whole sequence by the first point
if strcmp(opt.Normalize, 'Yes') == 1
    % Reshape the complex_data matrix into a 2D array for vectorized processing
    reshapedData = reshape(complex_data, [], size(complex_data, 4));
    % Compute the first phase for each signal group
    firstPhase = angle(reshapedData(:, 1));
    % Compute the manipulated signals
    manipulatedSignals = abs(reshapedData) .* exp(1i * (angle(reshapedData) - firstPhase));
    % Reshape the manipulated signals back into the original matrix shape
    complex_data = reshape(manipulatedSignals, matrixSize);
end
%%
% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)
complex_data_reoriented = write_volume(complex_data, input(1).nifti_header, 'Axial');
info = niftiinfo(files_in.In1{1});


nifti_header_output = info;
nifti_header_output.Filename = files_out.In1{1};
nifti_header_output.Filemoddate = char(datetime('now'));
[complex_data_reoriented, FinalMat] = CropNifti(complex_data_reoriented, nifti_header_output.Transform.T');
nifti_header_output.Datatype = class(complex_data_reoriented);
nifti_header_output.Transform = affine3d(FinalMat');
nifti_header_output.ImageSize = size(complex_data_reoriented); 
nifti_header_output.PixelDimensions = info.PixelDimensions(1:length(nifti_header_output.ImageSize));
nifti_header_output.MultiplicativeScaling = 1;

% % save the new .nii file
 niftiwrite(complex_data_reoriented, files_out.In1{1}, nifti_header_output);
 
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
