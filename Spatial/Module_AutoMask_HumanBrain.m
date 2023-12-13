function [files_in,files_out,opt] = Module_AutoMask(files_in,files_out,opt)

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
    module_option(:,4)   = {'OutOfTheMask','NaN'};
    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
    module_option(:,9)   = {'output_filename_ext','_AutoMasked'};
    module_option(:,10)   = {'output_ROIname','brain'};
    module_option(:,11)  = {'Output_orientation','Scan'};
    module_option(:,12)  = {'Crop_output','Yes'};
    module_option(:,13)  = {'threshold','55'};
    
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
    user_parameter(:,1)   = {'Description','Text','','','', '','Automatically skull-stripped MRI and MRF acquisitions, using thresholding.'}  ;
    user_parameter(:,2)   = {'Select the scan to mask','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'   .Output filename extension','char','','output_filename_ext','','',...
        {'Specify the string to be added to the first filename.'
        'Default filename extension is ''_AutoMasked''.'}'};
    user_parameter(:,4)   = {'   .Output ROI name','char','','output_ROIname','','',...
        {'Specify the saved ROI name.'
        'Default ROI name is ''brain''.'}'};
    user_parameter(:,5)   = {'   .Output orientation','cell',{'Scan', 'ROI'},'Output_orientation','','',...
        {'Specify the output orientation'
        '--> Output orienation = Scan'
        '--> Output orientation = ROI'
        }'};
    user_parameter(:,6) = {'   .Value out of the mask','char','','OutOfTheMask','','',''};
    user_parameter(:,7) = {'   .Crop empty slices','cell',{'Yes', 'No'},'Crop_output','','',...
        {'User can decide to crop empty slices'
        'for example, if the ROI selected is much smaller that the scan, the output masked scan can have (or not) the same number of slice than the ROI'
        'The default value is ''Yes''' 
        }'};
     user_parameter(:,8) = {'   .Thresholding (graylevel)','char','','threshold','','','threshold value used for cropping skull.'};
   
    
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
    Masked.Group = opt.Table_out.Group(1);
    Masked.Patient = opt.Table_out.Patient(1);
    Masked.Tp = opt.Table_out.Tp(1);
    Masked.Path = categorical(cellstr([opt.folder_out, filesep]));
    seqname = categorical(cellstr([char(opt.Table_out.SequenceName), opt.output_filename_ext]));
    Masked.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(seqname)]));
    Masked.Type = opt.Table_out.Type(1);
    Masked.IsRaw = categorical(0);   
    Masked.SequenceName = seqname;
    
    opt.Table_out = struct2table(Masked);
    files_out.In1{1} = [char(opt.Table_out.Path(1)), char(opt.Table_out.Filename(1)) '.nii'];
    
    ROIstruc.Group = opt.Table_out.Group(1);
    ROIstruc.Patient = opt.Table_out.Patient(1);
    ROIstruc.Tp = opt.Table_out.Tp(1);
    ROIstruc.Path = categorical(cellstr([opt.folder_out, filesep]));
    seqname = categorical(cellstr(opt.output_ROIname));
    ROIstruc.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(seqname)]));
    ROIstruc.Type = categorical(cellstr('ROI'));
    ROIstruc.IsRaw = categorical(0);  
    ROIstruc.SequenceName = seqname;

    opt.Table_out = [opt.Table_out; struct2table(ROIstruc)];
    files_out.In1{2} = [char(opt.Table_out.Path(2)), char(opt.Table_out.Filename(2)) '.nii'];
end





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Mask:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
if isfield(files_in, 'In2')
    input(2).nifti_header = spm_vol(files_in.In2{1});
end


if strcmp(opt.Output_orientation, 'Scan') %|| ~isfield(files_in, 'ROI')
    ref_scan = 1;
else
    ref_scan = 2;
end

Scan = read_volume(input(1).nifti_header, input(ref_scan).nifti_header, 0, 'Axial');
Types{1} = class(Scan);

[rows, columns, numberOfSlices, nTR] = size(Scan);

% Initialize variables to store results
brainROI_stack = zeros(rows, columns, numberOfSlices);
skullStrippedImage_stack = zeros(rows, columns, numberOfSlices);

% Loop through each slice
for sliceIdx = 1:numberOfSlices
    % Extract current slice
    grayImage_slice = mat2gray(Scan(:,:,sliceIdx)) * 256;

    % Thresholding to Create Binary Image
    thresholdValue = str2double(opt.threshold);
    binaryImage_slice = grayImage_slice > thresholdValue;
    binaryImage_slice = imclearborder(binaryImage_slice);

    % Extract Largest Blobs (Skull and Brain)
    binaryImage_slice = bwareafilt(binaryImage_slice, 2,  'largest');
    binaryImage_slice = imopen(binaryImage_slice, true(5));
    binaryImage_slice = bwareafilt(binaryImage_slice, 1,  'largest');
    binaryImage_slice = imfill(binaryImage_slice, 'holes');
    binaryImage_slice = imdilate(binaryImage_slice, true(5));

    % Mask Out Skull from Original Grayscale Image
    skullFreeImage_slice = grayImage_slice;
    skullFreeImage_slice(~binaryImage_slice) = 0;

    % Store results for the current slice
    brainROI_stack(:,:,sliceIdx) = binaryImage_slice;
    skullStrippedImage_stack(:,:,sliceIdx) = skullFreeImage_slice;
end

% % Display Brain ROI on Multiple Slices
% figure;
% montage(brainROI_stack, 'DisplayRange', [0 1]);
% title('Brain ROI - Multiple Slices');
% 
% % Display Skull Stripped Image on Multiple Slices
% figure;
% montage(skullStrippedImage_stack, 'DisplayRange', [0 255]);
% title('Skull Stripped Image - Multiple Slices');
% 

ROI = brainROI_stack;
CommonType = FindCommonDatatype(Types);       
Scan = cast(Scan, CommonType);
ROI = cast(ROI, CommonType);

if length(size(Scan)) > length(size(ROI))
    OutputImages = Scan .* repmat(ROI, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]); 
    ROItmp = repmat(ROI, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]);
    OutputImages(~ROItmp) = str2double(opt.OutOfTheMask);
else
    OutputImages =  repmat(Scan, [1 1 1 size(Scan,4) size(Scan,5) size(Scan,6) size(Scan,7)]) .* ROI;   
    OutputImages(~ROI) = str2double(opt.OutOfTheMask);

end

% OutputImages = skullStrippedImage_stack;


% transform the OutputImages matrix in order to match to the nii header of the
% first input (rotation/translation)
if ~exist('OutputImages_reoriented', 'var')
    OutputImages_reoriented = write_volume(OutputImages, input(ref_scan).nifti_header, 'Axial');
end

info = niftiinfo(files_in.(['In' num2str(ref_scan)]){1});

nifti_header_output = info;
nifti_header_output.Filename = files_out.In1{1};
nifti_header_output.Filemoddate = char(datetime('now'));

% since the option 'Crop_ouput' has been added afterwards, we first have to
% check if the option exist in order to exectute old pipeline.
if ~isfield(opt, 'Crop_output') || isfield(opt, 'Crop_output') && strcmp(opt.Crop_output, 'Yes')
    [OutputImages_reoriented, FinalMat] = CropNifti(OutputImages_reoriented, nifti_header_output.Transform.T');
    nifti_header_output.Transform = affine3d(FinalMat');
 
end

nifti_header_output.Datatype = class(OutputImages_reoriented);
nifti_header_output.ImageSize = size(OutputImages_reoriented); 
nifti_header_output.PixelDimensions = info.PixelDimensions(1:length(nifti_header_output.ImageSize));
nifti_header_output.MultiplicativeScaling = 1;

% % save the new .nii file
niftiwrite(OutputImages_reoriented, files_out.In1{1}, nifti_header_output);


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

%% ROI SAVING
ROI_reoriented = write_volume(ROI, input(ref_scan).nifti_header, 'Axial');

info2 = niftiinfo(files_in.In1{1});

info2.Filename = files_out.In1{2};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(ROI);
info2.ImageSize = info2.ImageSize(1:3);
info2.PixelDimensions = info2.PixelDimensions(1:3);


niftiwrite(ROI_reoriented, files_out.In1{2}, info2)

