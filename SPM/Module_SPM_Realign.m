function [files_in,files_out,opt] = Module_SPM_Realign(files_in,files_out,opt)

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
    module_option(:,3)   = {'OutputSequenceName','Prefix'};
    module_option(:,4)   = {'Est_Quality','0.9'};
    module_option(:,5)   = {'Est_Separation','4'};
    module_option(:,6)   = {'Est_Smoothing','5'};
    module_option(:,7)   = {'Est_Num_Passes','Register to first'};
    module_option(:,8)   = {'Est_Interpolation','2nd Degree B-Spline'};
    module_option(:,9)   = {'Est_Wrapping','No wrap'};
    
    module_option(:,10)   = {'Res_Interpolation','4th Degree B-Spline'};
    module_option(:,11)   = {'Res_Wrapping','No wrap'};
    module_option(:,12)   = {'Res_Masking','Dont mask images'};
    module_option(:,13)   = {'Res_output_filename_prefix','r'};
    module_option(:,14)   = {'RefInput',1};
    module_option(:,15)   = {'InputToReshape',1};
    module_option(:,16)   = {'Table_in', table()};
    module_option(:,17)   = {'Table_out', table()};
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
        {
        'This routine realigns a time-series of images acquired from the same subject using a least squares approach and a 6 parameter (rigid body) spatial transformation/* \cite{friston95a}*/.'
        'The first image is used as a reference to which all subsequent images are realigned.'
        ''
        'The aim is primarily to remove movement artefact in fMRI and PET time-series (or more generally longitudinal studies).' 
    ''
    }'};
    user_parameter(:,2)   = {'Image to realign','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please selelct the 4D scan to realign'};
    user_parameter(:,3)   = {'Parameters','','','','','',''};
    user_parameter(:,4)   = {'    Estimation Options','', '','','','',...
        'Various registration options. If in doubt, simply keep the default values.'};
    user_parameter(:,5)   = {'       .Quality','numeric','','Est_Quality','','',...
        'Quality versus speed trade-off.  Highest quality (1) gives most precise results, whereas lower qualities gives faster realignment. The idea is that some voxels contribute little to the estimation of the realignment parameters. This parameter is involved in selecting the number of voxels that are used.'};
    user_parameter(:,6)  = {'       .Separation','numeric','','Est_Separation','','',...
        'The separation (in mm) between the points sampled in the reference image.  Smaller sampling distances gives more accurate results, but will be slower.'};
    user_parameter(:,7)  = {'       .Histogram Smoothing','numeric','','Est_Smoothing','','',...
        {'The FWHM of the Gaussian smoothing kernel (mm) applied to the images before estimating the realignment parameters.'
        'Human MRI images typically use a 5 mm kernel.'
        'Rat MRI images typically use a 1 mm kernel.'}
        };
    user_parameter(:,8)  = {'       .Num Passes','cell',{'Register to first','Register to mean'},'Est_Num_Passes','','',...
        {'Register to first: Images are registered to the first image in the series. '
         'Register to mean:   A two pass procedure is used in order to register the images to the mean of the images after the first realignment.'
         ''
        'MRI images are typically registered to the first image.  The more accurate way would be to use a two pass procedure, but this probably wouldn''t improve the results so much and would take twice as long to run.'
        'PET images are typically registered to the mean. This is because PET data are more noisy than fMRI and there are fewer of them, so time is less of an issue.'}
        };
    user_parameter(:,9)  = {'       .Interpolation','cell',{'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'},'Est_Interpolation','','',...
        {'The method by which the images are sampled when estimating the optimum transformation.'
        'Higher degree interpolation methods provide the better interpolation, but they are slower because they use more neighbouring voxels /* \cite{thevenaz00a,unser93a,unser93b}*/. '}
        };
    user_parameter(:,10)  = {'       .Wrapping','cell',{'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'},'Est_Wrapping','','',...
        {'This indicates which directions in the volumes the values should wrap around in.  For example, in MRI scans, the images wrap around in the phase encode direction, so (e.g.) the subject''s nose may poke into the back of the subject''s head. These are typically:'
        '    No wrapping - for PET or images that have already been spatially transformed. Also the recommended option if you are not really sure.'
        '    Wrap in  Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'}
        }; 
    user_parameter(:,11)  = {'    Reslice options','','','','','',...
        'These images are resliced to the same dimensions, voxel sizes, orientation etc as the space defining image.'};
    user_parameter(:,12)  = {'       .Interpolation','cell',{'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'},'Res_Interpolation','','',...
        'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation is slower because it uses more neighbours.'};
    user_parameter(:,13)  = {'       .Wrapping','cell',{'No wrap','Wrap X', 'Wrap Y', 'Wrap X&Y', 'Wrap Z', 'Wrap X&Z', 'Wrap Y&Z', 'Wrap X,Y&Z'},'Res_Wrapping','','',...
        'These are typically: No wrapping - for PET or images that have already been spatially transformed.  Wrap in  Y  - for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'};
    user_parameter(:,14)  = {'       .Masking','cell',{'Mask images', 'Dont mask images'},'Res_Masking','','',...
       'Because of subject motion, different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).'};         
    user_parameter(:,15)  = {'       .Filename prefix','char','','Res_output_filename_prefix','','',...
        'Specify the string to be prepended to the filenames of the resliced image file(s). Default prefix is ''r''.'};

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


if strcmp(files_out, '')
    [Path_In1, Name_In1, ~] = fileparts(files_in.In1{1});
    tags1 = opt.Table_in(opt.Table_in.Path == [Path_In1, filesep],:);
    tags1 = tags1(tags1.Filename == Name_In1,:);
    assert(size(tags1, 1) == 1);
    tags_out_In1 = tags1;
    tags_out_In1.IsRaw = categorical(0);
    tags_out_In1.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_In1.SequenceName = categorical(cellstr([opt.Res_output_filename_prefix, char(tags_out_In1.SequenceName)]));
    tags_out_In1.Filename = categorical(cellstr([char(tags_out_In1.Patient), '_', char(tags_out_In1.Tp), '_', char(tags_out_In1.SequenceName)]));
    f_out = [char(tags_out_In1.Path), char(tags_out_In1.Patient), '_', char(tags_out_In1.Tp), '_', char(tags_out_In1.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    opt.Table_out = tags_out_In1;
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_Coreg_Est:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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

info = niftiinfo(files_in.In1{1});

if length(info.ImageSize) ~= 4
    error(['The file ', files_in.In1{1}, ' is not a 4D volume.'])
end

% First duplicate the source scan using the prefix string (user-defined)
% Otherwise the SPM_realign mayl overwrite the file!!
[~, In1_name, ~] = fileparts(files_in.In1{1});
copyfile(files_in.In1{1},  fullfile(char(opt.Table_out.Path), [In1_name, '.nii']));

header = spm_vol(fullfile(char(opt.Table_out.Path), [In1_name, '.nii']));
if numel(header) > 1
     data = cell([numel(header),1]);
    for j=1:numel(header)
        data{j}= [fullfile(char(opt.Table_out.Path), [In1_name, '.nii']), ',', num2str(j)];
    end
end


matlabbatch{1}.spm.spatial.realign.estwrite.data = {data};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = str2double(opt.Est_Quality);
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = str2double(opt.Est_Separation);
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = str2double(opt.Est_Smoothing);

switch opt.Est_Num_Passes
    case 'Register to first'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    case 'Register to mean'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
end

% Type of interpolation
switch opt.Est_Interpolation
    case 'Trilinear'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 7;
end

%  Type of Warpping
switch opt.Est_Wrapping
    case 'No wrap'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    case 'Warp X'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [1 1 1];
end

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';


% in MP3 do not used the mean scan (if generated). There the roptions.which
% is always set to [2 0]
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 0];

switch opt.Res_Interpolation
    case 'Nearest neighbour'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 0;
    case 'Trilinear'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 7;
end

switch opt.Res_Wrapping
    case 'No wrap'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    case 'Warp X'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [1 0 0];
    case 'Warp Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 1 0];
    case 'Warp X&Y'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [1 1 0];
    case 'Warp Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 1];
    case 'Warp X&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [1 0 1];
    case 'Warp Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 1 1];
    case 'Warp X,Y&Z'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [1 1 1];
end

switch opt.Res_Masking
    case 'Dont mask images'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 0;
    case 'Mask images'
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;

end

matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = opt.Res_output_filename_prefix;

jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});


%rename realigned file
movefile(fullfile(char(opt.Table_out.Path), [matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix, In1_name, '.nii']), files_out.In1{1});
%generate the jsonfile
jsonfile = strrep(files_out.In1{1}, '.nii', '.json');
copyfile(strrep(files_in.In1{1},'.nii','.json'),  jsonfile)
%% Json Processing
J = ReadJson(jsonfile);
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 
WriteJson(J, jsonfile)


% delete tmp files
delete(fullfile(char(opt.Table_out.Path), [In1_name, '.nii']));
delete(fullfile(char(opt.Table_out.Path), [In1_name, '.mat']));
delete(fullfile(char(opt.Table_out.Path), ['rp_', In1_name, '.txt']));

clear matlabbatch 

    









