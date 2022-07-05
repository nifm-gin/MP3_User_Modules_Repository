function [files_in,files_out,opt] = Module_pipeline_to_generate_FA_and_MD_maps(files_in,files_out,opt)


%%
% questions restantes : 
%     noms des scans d'entrées
%     vérifier la doc avec Michel

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
    module_option(:,3)   = {'output_filename_ext_MD','MD'};
    module_option(:,4)   = {'output_filename_ext_FA','FA'};
    module_option(:,5)   = {'OutputSequenceName','AllName'};
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
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {
        'Extraction of Fractional anisotropy (FA) and Mean diffusivity (MD) maps Diffusion-weighted images (DWI)'
        'Using MRtrix we can extract 2 diffusion tensor parameters : Fractional anisotropy (FA) and Mean diffusivity (MD)'
        ''
        'This module needs One DWI as input and eventually a second one acquired in opposite direction '
        'For instance input1 is acquired in anteroposterior and input2 in posteroanterior'
        'Output --> 2 parametrics maps: FA and MD '
        'This module correpond to a custom pipeline'
        '    designed by Arnaud Attyé, MD, CHU-Grenoble Alpes, France'
        '    Coded by Veronica Munoz Ramirez and adapted to MP3 by Benjamin Lemasson, Grenoble Institut of Neurosciences, France'
        ''
        'Briefly, this pipeline is composed of :'
        '     --> A DWI denoising (MRtrix : dwidenoise)'
        '     --> A Gibbs ringing artifacts removal (MRtrix : mrdegibbs)'
        '     --> A DWI distorsion correction (MRtrix : dwifslpreproc)'
        '     --> A Mask estimation (MRtrix : dwi2mask)'
        '     --> A DWI upsampling (MRtrix : mrgrid)'
        '     --> A Tensor-derived parameter maps generation (MRtrix : dwi2tensor & tensor2metric)'
        ''
        'WARNING'
        '   MRtrix3.0 needs to be installed (https://www.mrtrix.org/) in /usr/local/mrtrix3'
        '   FSL 6.0.1 needs to be installed by the user first in /usr/local/fsl'
        }'};
    
    user_parameter(:,2)   = {'Select the Raw DWI-scan','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Select the DWI-scan acquired in opposite direction (Optional)','1Scan','','',{'SequenceName'}, 'Optional',''};
    
    user_parameter(:,4)   = {'Parameters','','','','', '', ''};
    user_parameter(:,5)   = {'   .Output filename','char','FA','output_filename_ext_FA','','',...
        {'Specify the name of the fist output scan : the Fractional anisotropy.'
        'Default filename is ''FA''.'}'};
    user_parameter(:,6)   = {'   .Output filename','char','MD','output_filename_ext_MD','','',...
        {'Specify the name of the second output scan : the Mean diffusivity.'
        'Default filename is ''MD''.'}'};
    
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

opt.NameOutFiles = {opt.output_filename_ext_FA, opt.output_filename_ext_MD};

if isempty(files_out)
    for i=1:length(opt.NameOutFiles)
        table_out_tmp = opt.Table_in(1,:);
        table_out_tmp.IsRaw = categorical(0);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        if strcmp(opt.OutputSequenceName, 'AllName')
            table_out_tmp.SequenceName = categorical(cellstr(opt.NameOutFiles{i}));
        elseif strcmp(opt.OutputSequenceName, 'Extension')
            table_out_tmp.SequenceName = categorical(cellstr([char(table_out_tmp.SequenceName), opt.NameOutFiles{i}]));
        end
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        opt.Table_out = [opt.Table_out; table_out_tmp];
    end
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


mrtrix_path = '/usr/local/mrtrix3/bin/';

setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('PATH', [getenv('PATH') ':/usr/local/mrtrix3/bin']);

setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be



% Since MRtrix needs bvec and bval files to generate FA and MD files we
% have to generate such files using .json files
%% first we will generate the bvecs and bvals for input1 and copy input1 in the tmp folder
DWI_json_data = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));
[~,NAME,~] = fileparts(files_in.In1{1});
%bvals
DWI_bvals = strrep(DWI_json_data.bvals.value, ' ', ',');
format long
if ischar(DWI_bvals)
    DWI_bvals = str2num(DWI_bvals(3:end-2)); %#ok<ST2NM>
else
    DWI_bvals = str2num(DWI_bvals{1,1}); %#ok<ST2NM>
end
    
fileID = fopen(fullfile(opt.folder_out, [NAME '.bvals']),'w');
fprintf(fileID,'%5f ',DWI_bvals);
fclose(fileID);

%bvecs
DWI_bvecs = strrep(DWI_json_data.bvecs.value, ' ', ',');
DWI_bvecs = str2num(DWI_bvecs{1,1}); %#ok<ST2NM>
fileID = fopen(fullfile(opt.folder_out, [NAME '.bvecs']),'w');
fprintf(fileID, [repmat('%.16f ', 1, size(DWI_bvecs,2)-1), '%.20f\n'], DWI_bvecs');
fclose(fileID);

%copy input1 in the tmp folder
copyfile(files_in.In1{1},  fullfile(opt.folder_out, [NAME '.nii']))
DWI_file = fullfile(opt.folder_out, [NAME '.nii']);



% MR convert nifti to mif
disp('-------- MR convert nifti to mif ---------');
status = system([mrtrix_path 'mrconvert ' DWI_file ' ' fullfile(opt.folder_out, 'DWI.mif') ' -fslgrad ' strrep(DWI_file, '.nii', '.bvecs') ' ' strrep(DWI_file, '.nii', '.bvals')]);
if status ~= 0
    warning('FAIL - mrconvert input 1');
end
% DWI denoising
disp('-------- DWI denoising ---------');
status = system([mrtrix_path 'dwidenoise ' fullfile(opt.folder_out, 'DWI.mif ') fullfile(opt.folder_out, 'DWI_denoised.mif ') ' -noise ' fullfile(opt.folder_out, 'noise_DWI_opp.mif')]);
if status ~= 0
    warning('FAIL - dwidenoise input 1');
end

% Gibbs ringing artifacts removal
disp('-------- Gibbs ringing artifacts removal ---------');
status = system([mrtrix_path 'mrdegibbs ' fullfile(opt.folder_out, 'DWI_denoised.mif ') fullfile(opt.folder_out, 'DWI_Gibbs.mif ')]);
if status ~= 0
    warning('FAIL - mrdegibbs input 1');
end


%% Now we will generate the bvecs and bvals for input2 and processing (in present) : the DWI_opp_dir scan
if isfield(files_in, 'In2')
    DWI_opp_dir_json_data = spm_jsonread(strrep(files_in.In2{1}, '.nii', '.json'));
    [~,NAME,~] = fileparts(files_in.In2{1});
    %bvals
    DWI_opp_dir_bvals = strrep(DWI_opp_dir_json_data.bvals.value, ' ', ',');
    
    if ischar(DWI_opp_dir_bvals)
        DWI_opp_dir_bvals = str2num(DWI_opp_dir_bvals(3:end-2)); %#ok<ST2NM>
    else
        DWI_opp_dir_bvals = str2num(DWI_opp_dir_bvals{1,1}); %#ok<ST2NM>
    end
    
    fileID = fopen(fullfile(opt.folder_out, [NAME '.bvals']),'w');
    fprintf(fileID,'%5f ',DWI_opp_dir_bvals);
    fclose(fileID);
    
    %bvecs
    DWI_opp_dir_bvecs = strrep(DWI_opp_dir_json_data.bvecs.value, ' ', ',');
    DWI_opp_dir_bvecs = str2num(DWI_opp_dir_bvecs{1,1}); %#ok<ST2NM>
    fileID = fopen(fullfile(opt.folder_out, [NAME '.bvecs']),'w');
    fprintf(fileID, [repmat('%16f ', 1, size(DWI_opp_dir_bvecs,2)-1), '%20f\n'], DWI_opp_dir_bvecs');
    fclose(fileID);
    
    %copy input2 in the tmp folder
    copyfile(files_in.In1{1},  fullfile(opt.folder_out, [NAME '.nii']))
    DWI_opp_dir_file = fullfile(opt.folder_out, [NAME '.nii']);
    
    % input 1 processing : the DWI scan
    % MR convert nifti to mif
    disp('-------- MR convert nifti to mif ---------');
    status = system([mrtrix_path 'mrconvert ' DWI_file ' ' fullfile(opt.folder_out, 'DWI_opp_dir.mif') ' -fslgrad ' strrep(DWI_file, '.nii', '.bvecs') ' ' strrep(DWI_file, '.nii', '.bvals')]);
    if status ~= 0
        warning('FAIL - mrconvert APA input 2');
    end
    
    % DWI denoising
    disp('-------- DWI denoising ---------');
    status = system([mrtrix_path 'dwidenoise ' fullfile(opt.folder_out, 'DWI_opp_dir.mif ') fullfile(opt.folder_out, 'DWI_opp_dir_denoised.mif ') ' -noise ' fullfile(opt.folder_out, 'noise_DWI_opp_dir.mif')]);
    if status ~= 0
        warning('FAIL - dwidenoise APP input 2');
    end
    
    % Gibbs ringing artifacts removal
    disp('-------- Gibbs ringing artifacts removal ---------');
    status = system([mrtrix_path 'mrdegibbs ' fullfile(opt.folder_out, 'DWI_opp_dir_denoised.mif ') fullfile(opt.folder_out, 'DWI_opp_dir_Gibbs.mif ')]);
    if status ~= 0
        warning('FAIL - mrdegibbs APP input 2');
    end
    
    % DWI distorsion correction
    disp('-------- DWI distorsion correction ---------');
    system([mrtrix_path 'mrcat ' fullfile(opt.folder_out, 'DWI_Gibbs.mif ') fullfile(opt.folder_out, 'DWI_opp_dir_Gibbs.mif ') fullfile(opt.folder_out, 'all_DWIs.mif ') -axis 3']);
    if status ~= 0
        warning('FAIL - mrcat');
    end
    
    
    if exist([mrtrix_path 'dwipreproc'])~=0
        status = system([mrtrix_path 'dwipreproc all_DWIs.mif dwi_preproc.mif -pe_dir PA -rpe_all -eddy_options="--slm=linear" ']);
    elseif exist([mrtrix_path 'dwifslpreproc'])~=0
        tic
        status = system([mrtrix_path 'dwifslpreproc ' fullfile(opt.folder_out, 'all_DWIs.mif ') fullfile(opt.folder_out, 'dwi_preproc.mif ') '-pe_dir PA -rpe_all -eddy_options="--slm=linear" ']);
        toc
    else
        status = 0;
    end
    if status ~= 0
        warning('FAIL - dwipreproc');
    end
else
    movefile(fullfile(opt.folder_out, 'DWI_Gibbs.mif'), fullfile(opt.folder_out, 'dwi_preproc.mif'))
end



% Mask estimation
disp('-------- Mask estimation ---------');
status = system([mrtrix_path 'dwi2mask ' fullfile(opt.folder_out, 'dwi_preproc.mif ') fullfile(opt.folder_out, 'mask_orig.mif ')]);
if status ~= 0
    warning('FAIL - dwi2mask');
end

% DWI upsampling
disp('-------- DWI bias correction and upsampling ---------');
status = system([mrtrix_path 'mrgrid ' fullfile(opt.folder_out, 'dwi_preproc.mif ') 'regrid -vox 1 ' fullfile(opt.folder_out, 'dwi_preproc_upsampled.mif')]);

if status ~= 0
    warning('FAIL - upsampling');
end
status = system([mrtrix_path 'dwi2mask ' fullfile(opt.folder_out, 'dwi_preproc_upsampled.mif ')  fullfile(opt.folder_out, 'mask_upsampled.mif')]);
if status ~= 0
    warning('FAIL - mask upsampling');
end

% Tensor-derived parameter maps generation
disp('-------- Tensor-derived parameter maps generation ---------');
status = system([mrtrix_path 'dwi2tensor -mask ' fullfile(opt.folder_out, 'mask_orig.mif ') fullfile(opt.folder_out, 'dwi_preproc.mif ')  fullfile(opt.folder_out, 'dwi_tensor.mif')]);
if status ~= 0
    warning('FAIL - dwi2tensor');
end
status = system([mrtrix_path 'tensor2metric -adc ' files_out.In1{2} ' -fa ' files_out.In1{1} ' -mask ' fullfile(opt.folder_out, 'mask_orig.mif ') fullfile(opt.folder_out, 'dwi_tensor.mif')]);
if status ~= 0
    warning('FAIL - tensor2metric');
end

%% Save file and create Jsonfile
Json_to_update = ReadJson(strrep(files_in.In1{1}, '.nii', '.json'));

for i=1:length(opt.NameOutFiles)
    Json_to_save= KeepModuleHistory(Json_to_update, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename); 

    [path, name, ~] = fileparts(files_out.In1{i});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(Json_to_save, jsonfile)
end

% delete tmp files
files_to_delete = dir(fullfile(opt.folder_out, '*.mif'));
for i=1:length(files_to_delete)
    delete(fullfile(opt.folder_out, files_to_delete(i).name));
end
files_to_delete = dir(fullfile(opt.folder_out, '*.bvecs'));
for i=1:length(files_to_delete)
    delete(fullfile(opt.folder_out, files_to_delete(i).name));
end
files_to_delete = dir(fullfile(opt.folder_out, '*.bvals'));
for i=1:length(files_to_delete)
    delete(fullfile(opt.folder_out, files_to_delete(i).name));
end
if isfield(files_in, 'In2')
    delete(DWI_opp_dir_file)
end
delete(DWI_file)
