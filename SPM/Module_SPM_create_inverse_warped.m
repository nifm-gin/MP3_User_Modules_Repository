function [files_in,files_out,opt] = Module_SPM_create_inverse_warped(files_in,files_out,opt)

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
    module_option(:,3)   = {'Time_Steps','64'};
    module_option(:,4)   = {'Interpolation','Trilinear'};

    module_option(:,5)   = {'RefInput',1};
    module_option(:,6)   = {'InputToReshape',1};
    module_option(:,7)   = {'Table_in', table()};
    module_option(:,8)   = {'Table_out', table()};
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
    user_parameter(:,1)   = {'Description','Text','','','', '',...
        {'Create Inverse Warped as described in spm12'
        'Create inverse normalised versions of some image(s). The image that is inverse-normalised should be in alignment with the template (generated during the warping procedure).' 
        'Note that the results have the same dimensions as the ``flow fields'''', but are mapped to the original images via the affine transformations in their headers.'
        }}  ;
    user_parameter(:,2)   = {'Images','XScanOrXROI','','',{'SequenceName'}, 'Mandatory',{'Select the image(s) to be inverse normalised.'
        'These should be in alignment with the template image of the warping procedure (Run Dartel).'}};
    user_parameter(:,3)   = {'Flow fields','1Scan','','',{'SequenceName'}, 'Mandatory',{'The flow fields store the deformation information.'
        'The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'}};
   
    user_parameter(:,4)   = {'Parameters','','','','','',''};
    user_parameter(:,5)   = {'   .Time Steps','cell',{'1', '2','4','8','16','32','64','128','256','512'},'Time_Steps','','',...
        {'The number of time points used for solving the partial differential equations.'
        'Note that Jacobian determinants are not very accurate for very small numbers of time steps (less than about 16).'
        ''
        'One of the following options must be selected:'
        '1'
        '2'
        '4'
        '8'
        '16'
        '32'
        '64'
        '128'
        '256'
        '512'
        }'};
    user_parameter(:,6)   = {'   .Interpolation','cell',{'Nearest neighbour', 'Trilinear', '2nd Degree B-Spline', '3rd Degree B-Spline', '4th Degree B-Spline', '5th Degree B-Spline', '6th Degree B-Spline', '7th Degree B-Spline'},'Interpolation','','',...
        'The method by which the images are sampled when being written in a different space. Nearest Neighbour is fastest, but not normally recommended. It can be useful for re-orienting images while preserving the original intensities (e.g. an image consisting of labels). Trilinear Interpolation is OK for PET, or realigned and re-sliced fMRI. If subject movement (from an fMRI time series) is included in the transformations then it may be better to use a higher degree approach. Note that higher degree B-spline interpolation is slower because it uses more neighbours.'};
  
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

    opt.IsRaw = categorical(0);
    opt.Path = categorical(cellstr([opt.folder_out, filesep]));

    for i=1:length(files_in.In1)
        table_out_tmp = opt.Table_in(i,:);
        table_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        table_out_tmp.SequenceName = categorical(cellstr(['wu', char(table_out_tmp.SequenceName)]));
        table_out_tmp.Filename = categorical(cellstr([char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName)]));
        opt.Table_out = [opt.Table_out; table_out_tmp];
        
        f_out = [char(table_out_tmp.Path), char(table_out_tmp.Patient), '_', char(table_out_tmp.Tp), '_', char(table_out_tmp.SequenceName), '.nii'];
        files_out.In1{i} = f_out;
        
    end
end



%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_SPM_create_inverse_warped:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 0;
    case 'Trilinear'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 1;
    case '2nd Degree B-Spline'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 2;
    case '3rd Degree B-Spline'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 3;
    case '4th Degree B-Spline'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 4;
    case '5th Degree B-Spline'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 5;
    case '6th Degree B-Spline'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 6;
    case '7th Degree B-Spline'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 7;
end

switch opt.Time_Steps
    case  '1'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 0;
    case   '2'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 1;
    case   '4'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 2;
    case   '8'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 3;
    case    '16'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 4;
    case   '32'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 5;
    case   '64'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 6;
    case    '128'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 7;
    case    '256'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 8;
    case    '512'
        matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 9;
end

% First duplicate the source scan using the prefix string (user-defined)
% Otherwise the CoregEstimate will overwrite the file!!
tmp_flowfield_fullfilename =  [opt.Table_in.Properties.UserData.MP3_data_path, 'Tmp', filesep, char(opt.Table_in.Filename(end)), '.nii'];
% if this module is executed in a pipeline, is it possible that
% files_in.In2{1} is already in the Tmp --> if so, skip
if ~strcmp(files_in.In2{1}, tmp_flowfield_fullfilename)
    copyfile(files_in.In2{1}, tmp_flowfield_fullfilename)
end
matlabbatch{1}.spm.tools.dartel.crt_iwarped.flowfields =  {tmp_flowfield_fullfilename};

images_list = {};
for i=1:length(files_in.In1)
    if ~isempty(files_in.In1{i})
        %header = niftiinfo(files_in.In3{i});
        header = spm_vol(files_in.In1{i});
        % First duplicate all the 'Other' scans using the prefix string (user-defined)
        % Otherwise the CoregEstimate will overwrite the raw files!!

        
        tmp_images_fullfilename =  [opt.Table_in.Properties.UserData.MP3_data_path, 'Tmp', filesep, char(opt.Table_in.Filename(i)), '.nii'];
        % if this module is executed in a pipeline, is it possible that
        % files_in.In2{1} is already in the Tmp --> if so, skip
        if ~strcmp(files_in.In1{i},  tmp_images_fullfilename)
            copyfile(files_in.In1{i},  tmp_images_fullfilename)
        end
        if ~strcmp(strrep(files_in.In1{i},'.nii','.json'), strrep(tmp_images_fullfilename,'.nii','.json'))
            if exist(strrep(files_in.In1{i},'.nii','.json'), 'file')
                copyfile(strrep(files_in.In1{i},'.nii','.json'),  strrep(tmp_images_fullfilename,'.nii','.json'))
            end
        end
        % copyfile(files_in.In1{i},  files_out.In1{i})
        
        if numel(header) > 1
            for j=1:numel(header)
                images_list= [images_list, tmp_images_fullfilename];
            end
        else
            images_list = [images_list, tmp_images_fullfilename];
        end
    end
end
matlabbatch{1}.spm.tools.dartel.crt_iwarped.images = images_list';



% [SPMinter,SPMgraph,~] = spm('FnUIsetup','test',1);
jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});

% rename the output file in order to match with what user expect
for i=1:length(files_in.In1)
    output_fullfilename = [opt.Table_in.Properties.UserData.MP3_data_path, 'Tmp', filesep, char(opt.Table_out.Filename(i)), '.nii'];
    movefile(strcat(opt.Table_in.Properties.UserData.MP3_data_path, 'Tmp', filesep, ['w', char(opt.Table_in.Filename(i)),'_', char(opt.Table_in.Filename(end))],'.nii'), output_fullfilename);
    if exist(strrep(files_in.In1{i},'.nii','.json'), 'file')
        copyfile(strrep(files_in.In1{i},'.nii','.json'),  strrep(output_fullfilename,'.nii','.json') )
        %% Json processing
        J = ReadJson(strrep(output_fullfilename,'.nii','.json') );
        J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
        WriteJson(J, strrep(output_fullfilename,'.nii','.json') )
    end
end




