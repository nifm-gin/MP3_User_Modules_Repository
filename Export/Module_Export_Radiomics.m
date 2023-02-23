function [files_in,files_out,opt] = Module_Export_Radiomics(files_in,files_out,opt)

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
    module_option(:,7)   = {'Output_Folder','Radiomics'};
    module_option(:,8)   = {'NbBin',100};
    module_option(:,9)   = {'NbGrayLevels',100};
    module_option(:,10)   = {'Isotropic_Resampling_Val',1};
    module_option(:,11)   = {'AutomaticJobsCreation', 'No'};
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
         {'This module computes radiomic features from an image and a ROI. It is based on the radiomic features computation described in Vallières, M. et al. (2015). A radiomics model from joint FDG-PET and MRI texture features for the prediction of lung metastases in soft-tissue sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 5471-5496. doi:10.1088/0031-9155/60/14/5471 . (Code available at https://github.com/mvallieres/radiomics)'}};
    user_parameter(:,2)   = {'   .Scan','1Scan','','', {'SequenceName'},'Mandatory',...
         'Please select the scans that will be analysed'};
    user_parameter(:,3)   = {'   .ROI','1ROI','','',{'SequenceName'},'Mandatory',...
         'Please select the ROI that will be analysed'};
    user_parameter(:,4)   = {'   .Output folder','char','','Output_Folder','', '','the name of the folder where will be saved the metrics, inside your project folder.'};
    user_parameter(:,5)   = {'   .Number of bins for histogram computation ?','numeric','','NbBin','', '',''};
    user_parameter(:,6)   = {'   .Number of gray levels for image quantization ?','numeric','','NbGrayLevels','', '',''};
    user_parameter(:,7)   = {'   .Isotropic Resampling : Resolution of the image (in mm3) ?','numeric','','Isotropic_Resampling_Val','', '',''};
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

UPat = unique(opt.Table_in.Patient);
UTp = unique(opt.Table_in.Tp);

for i=1:length(UTp)
    DBTp = opt.Table_in(opt.Table_in.Tp == UTp(i),:);
    for j=1:length(UPat)
        DBPatTp = DBTp(DBTp.Patient == UPat(j),:);
        ROIDB = DBPatTp(DBPatTp.Type == 'ROI',:);       
        if isempty(ROIDB)
            continue
        elseif height(ROIDB)>1
            warning('Plus d''une ROI trouvée, pas de figure générée')
            continue
        end
        ROIFilename = [char(ROIDB.Path), char(ROIDB.Filename), '.nii'];
        ROI_h = spm_vol(ROIFilename);
        ROI_vol = read_volume(ROI_h, ROI_h, 0, 'Axial');
        
        DBScan = DBPatTp(DBPatTp.Type == 'Scan',:); 
        if isempty(DBScan)
                continue
        elseif height(DBScan)>1
                warning('Plus d''un scan trouvé, pas de figure générée pour ce scan')
                continue
        end
        ScanFilename = [char(DBScan.Path), char(DBScan.Filename), '.nii'];
        scan_h = spm_vol(ScanFilename);
        scan_vol = read_volume(scan_h, ROI_h, 0, 'Axial');
        
        hdr = spm_read_hdr(ROIFilename);
        pixdim_x_y = hdr.dime.pixdim(2);
        pixdim_z = hdr.dime.pixdim(4);
        %info = niftiinfo(ROIFilename)
        %pixdim_x_y = info.PixelDimensions(1);
        %pixdim_z = info.PixelDimensions(3);
        
        [ROIonly] = prepareVolume(scan_vol,ROI_vol,'CT',pixdim_x_y,pixdim_z,1,opt.Isotropic_Resampling_Val,'Global');
        [globalTextures] = getGlobalTextures(ROIonly,opt.NbBin);
        %[eccentricity] = getEccentricity(ROIonly,pixdim_x_y,pixdim_z);
        %[sizeROI] = getSize(ROIonly,pixdim_x_y,pixdim_z);
        %[solidity] = getSolidity(ROIonly,pixdim_x_y,pixdim_z);
        %[volume] = getVolume(ROIonly,pixdim_x_y,pixdim_z);
        
        [ROIonly,levels] = prepareVolume(scan_vol,ROI_vol,'CT',pixdim_x_y,pixdim_z,1,opt.Isotropic_Resampling_Val,'Matrix','Uniform',opt.NbGrayLevels);
        [GLCM] = getGLCM(ROIonly,levels); 
        [glcmTextures] = getGLCMtextures(GLCM);
        glcmTextures.Variance_glcm = glcmTextures.Variance;
        glcmTextures = rmfield(glcmTextures, 'Variance');
        [GLRLM] = getGLRLM(ROIonly,levels); 
        [glrlmTextures] = getGLRLMtextures(GLRLM);
        [GLSZM] = getGLSZM(ROIonly,levels); 
        [glszmTextures] = getGLSZMtextures(GLSZM);
        glszmTextures.GLN_glszm = glszmTextures.GLN;
        glszmTextures = rmfield(glszmTextures, 'GLN');
        glszmTextures.GLV_glszm = glszmTextures.GLV;
        glszmTextures = rmfield(glszmTextures, 'GLV');
        [NGTDM,countValid] = getNGTDM(ROIonly,levels); 
        [ngtdmTextures] = getNGTDMtextures(NGTDM,countValid);
        ngtdmTextures.Contrast_ngtdm = ngtdmTextures.Contrast;
        ngtdmTextures = rmfield(ngtdmTextures, 'Contrast');
        
        parameters = [fieldnames(globalTextures)', fieldnames(glcmTextures)', fieldnames(glrlmTextures)', fieldnames(glszmTextures)', fieldnames(ngtdmTextures)'];
        current_cvs_table = table;
        current_cvs_table(1,1:size([{'PatientName'}, {'TimePoint'}, {'Reference'},  {'Scan_name'},  {'ROI_name'},  parameters],2)) = ...
                    [cellstr(UPat(j)), cellstr(UTp(i)), cellstr(ROIDB.SequenceName),  cellstr(DBScan.SequenceName), cellstr(ROIDB.SequenceName), num2cell(NaN(1, size(parameters,2)))];
        parameters = clean_variable_name(parameters,1);
        current_cvs_table.Properties.VariableNames = [{'PatientName'}, {'TimePoint'}, {'Reference'},  {'Scan_name'},  {'ROI_name'},  parameters];
        %current_cvs_table(:,7:end) = num2cell(computed_values);
        
        for aaa=1:length(parameters)
            param = parameters(aaa);
            if isfield(globalTextures, param{1})
                current_cvs_table = setfield(current_cvs_table, param{1}, getfield(globalTextures, param{1}));
            elseif isfield(glcmTextures, param{1})
                current_cvs_table = setfield(current_cvs_table, param{1}, getfield(glcmTextures, param{1}));
            elseif isfield(glrlmTextures, param{1})
                current_cvs_table = setfield(current_cvs_table, param{1}, getfield(glrlmTextures, param{1}));
            elseif isfield(glszmTextures, param{1})
                current_cvs_table = setfield(current_cvs_table, param{1}, getfield(glszmTextures, param{1}));
            elseif isfield(ngtdmTextures, param{1})
                current_cvs_table = setfield(current_cvs_table, param{1}, getfield(ngtdmTextures, param{1}));
            end
        end
        
        if ~exist('cvs_table', 'var')
            cvs_table = current_cvs_table;
        else
            cvs_table = outerjoin(cvs_table,current_cvs_table,'MergeKeys', true);
        end
        
    end
end

folder_out_silesep = strsplit(opt.folder_out, filesep);
folder_out_silesep(end) = {'Data_to_export'};
folder_out_silesep(end+1) = {opt.Output_Folder};
folder_out_silesep(end+1) = {''};

if exist(strjoin(folder_out_silesep, filesep),'dir') ~= 7
    [status, ~, ~] = mkdir(strjoin(folder_out_silesep, filesep));
    if status == false
        error('Cannot create the folder to save the csv files.')
    end
end
f_out = [strjoin(folder_out_silesep, filesep), opt.Output_Folder, '_radiomicFeatures.csv'];


writetable(cvs_table, f_out)


