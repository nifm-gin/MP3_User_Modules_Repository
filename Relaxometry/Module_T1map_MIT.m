function [files_in,files_out,opt] = Module_T1map_MIT(files_in,files_out,opt)
% This is a template file for "brick" functions in NIAK.
%
% SYNTAX:
% [IN,OUT,OPT] = PSOM_TEMPLATE_BRICK(IN,OUT,OPT)
%
% _________________________________________________________________________
% INPUTS:
%
% IN
%   (string) a file name of a 3D+t fMRI dataset .
%
% OUT
%   (structure) with the following fields:
%       flag_test
%   CORRECTED_DATA
%       (string, default <BASE NAME FMRI>_c.<EXT>) File name for processed
%       data.
%       If OUT is an empty string, the name of the outputs will be
%       the same as the inputs, with a '_c' suffix added at the end.
%
%   MASK
%       (string, default <BASE NAME FMRI>_mask.<EXT>) File name for a mask
%       of the data. If OUT is an empty string, the name of the
%       outputs will be the same as the inputs, with a '_mask' suffix added
%       at the end.
%
% OPT
%   (structure) with the following fields.
%
%   TYPE_CORRECTION
%      (string, default 'mean_var') possible values :
%      'none' : no correction at all
%      'mean' : correction to zero mean.
%      'mean_var' : correction to zero mean and unit variance
%      'mean_var2' : same as 'mean_var' but slower, yet does not use as
%      much memory).
%
%   FOLDER_OUT
%      (string, default: path of IN) If present, all default outputs
%      will be created in the folder FOLDER_OUT. The folder needs to be
%      created beforehand.
%
%   FLAG_VERBOSE
%      (boolean, default 1) if the flag is 1, then the function prints
%      some infos during the processing.
%
%   FLAG_TEST
%      (boolean, default 0) if FLAG_TEST equals 1, the brick does not do
%      anything but update the default values in IN, OUT and OPT.
%
% _________________________________________________________________________
% OUTPUTS:
%
% IN, OUT, OPT: same as inputs but updated with default values.
%
% _________________________________________________________________________
% SEE ALSO:
% NIAK_CORRECT_MEAN_VAR
%
% _________________________________________________________________________
% COMMENTS:
%
% _________________________________________________________________________
% Copyright (c) <NAME>, <INSTITUTION>, <START DATE>-<END DATE>.
% Maintainer : <EMAIL ADDRESS>
% See licensing information in the code.
% Keywords : PSOM, documentation, template, brick

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization and syntax checks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the module's parameters with default values
if isempty(opt)
    
    % define every option needed to run this module
    %%   % define every option needed to run this module
    % --> module_option(1,:) = field names
    % --> module_option(2,:) = defaults values
    module_option(:,1)   = {'folder_out',''};
    module_option(:,2)   = {'flag_test',true};
    module_option(:,3)   = {'first_scan',1};
    module_option(:,4)   = {'M0_map','No'};
    module_option(:,5)   = {'IE_map','No'};
    module_option(:,6)   = {'T1_output_filename','T1Map'};
    module_option(:,7)   = {'M0_output_filename','M0map'};
    module_option(:,8)   = {'IE_output_filename','IEmap'};
    module_option(:,9)   = {'OutputSequenceName','AllName'};
    module_option(:,10)   = {'RefInput',1};
    module_option(:,11)   = {'InputToReshape',1};
    module_option(:,12)   = {'Table_in', table()};
    module_option(:,13)   = {'Table_out', table()};
    opt.Module_settings = psom_struct_defaults(struct(),module_option(1,:),module_option(2,:));
    
    % list of everything displayed to the user associated to their 'type'
    % --> user_parameter(1,:) = user_parameter_list
    % --> user_parameter(2,:) = user_parameter_type
    % --> user_parameter(3,:) = parameter_default
    % --> user_parameter(4,:) = psom_parameter_list
    % --> user_parameter(5,:) = Scans_Input_DOF (degree-of-freedom)
    % --> user_parameter(6,:) = Help : text data which describe the parameter (it
    % will be display to help the user)
    user_parameter(:,1)   = {'Description','Text','','','','',...
        {'Generate a T1map from a Multi-Inversion Time scan'}'};
    user_parameter(:,2)   = {'Select a MTI_T1 scan as input','1Scan','','',{'SequenceName'}, 'Mandatory',''};
    user_parameter(:,3)   = {'Parameters','','','','', '',''};
    user_parameter(:,4)   = {'   .Output filename','char','T1Map','T1_output_filename','','',...
        {'Specify the name of the T1map generated'
        'Default filename is ''T1map''.'}'};
    user_parameter(:,5)   = {'   .First Scan','numeric','1','first_scan','', '','Please, select the first dynamic to use; by default the all dynamic is used'};
    user_parameter(:,6)   = {'   .Save the M0 map ?','cell',{'Yes', 'No'},'M0_map','', '',''};
    user_parameter(:,7)  =  {'       .M0 Filename ','char','','M0_output_filename','','',...
        'Specify the name of the M0 map'};
    user_parameter(:,8)   = {'   .IE map ?','cell',{'Yes', 'No'},'IE_map','', '',''};
    user_parameter(:,9)  =  {'       .IE Filename ','char','','IE_output_filename','','',...
        'Specify the name of the IE map'};
    
    
    VariableNames = {'Names_Display', 'Type', 'Default', 'PSOM_Fields', 'Scans_Input_DOF', 'IsInputMandatoryOrOptional','Help'};
    opt.table = table(user_parameter(1,:)', user_parameter(2,:)', user_parameter(3,:)', user_parameter(4,:)', user_parameter(5,:)', user_parameter(6,:)',user_parameter(7,:)', 'VariableNames', VariableNames);
    
    % So far no input file is selected and therefore no output
    %
    % The output file will be generated automatically when the input file
    % will be selected by the user
    opt.NameOutFiles = {'ASL_InvEff'};
    
    files_in.In1 = {''};
    files_out.In1 = {''};
    return
    
end
%%%%%%%%
% debug
if isfield(opt, 'output_filename_ext')
    opt.output_filename = opt.output_filename_ext;
end

opt.NameOutFiles = {opt.T1_output_filename};


if isempty(files_out)
    opt.Table_out = opt.Table_in;
    opt.Table_out.IsRaw = categorical(0);
    opt.Table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
    if strcmp(opt.OutputSequenceName, 'AllName')
        opt.Table_out.SequenceName = categorical(cellstr(opt.T1_output_filename));
    elseif strcmp(opt.OutputSequenceName, 'Extension')
        opt.Table_out.SequenceName = categorical(cellstr([char(opt.Table_out.SequenceName), opt.T1_output_filename]));
    end
    opt.Table_out.Filename = categorical(cellstr([char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName)]));
    f_out = [char(opt.Table_out.Path), char(opt.Table_out.Patient), '_', char(opt.Table_out.Tp), '_', char(opt.Table_out.SequenceName), '.nii'];
    files_out.In1{1} = f_out;
    
    if strcmp(opt.M0_map, 'Yes')
        tmp = opt.Table_out(1,:);
        tmp.SequenceName = categorical(cellstr(opt.M0_output_filename));
        tmp.Filename = categorical(cellstr([char(tmp.Patient), '_', char(tmp.Tp), '_', char(tmp.SequenceName)]));
        f_out = [char(tmp.Path), char(tmp.Patient), '_', char(tmp.Tp), '_', char(tmp.SequenceName), '.nii'];
        files_out.In2{1} = f_out;
        opt.Table_out = [opt.Table_out; tmp];
    end
    
    if strcmp(opt.IE_map, 'Yes')
        tmp = opt.Table_out(1,:);
        tmp.SequenceName = categorical(cellstr(opt.IE_output_filename));
        tmp.Filename = categorical(cellstr([char(tmp.Patient), '_', char(tmp.Tp), '_', char(tmp.SequenceName)]));
        f_out = [char(tmp.Path), char(tmp.Patient), '_', char(tmp.Tp), '_', char(tmp.SequenceName), '.nii'];
        files_out.In3{1} = f_out;
        opt.Table_out = [opt.Table_out; tmp];
    end
    
end

%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('T1map_MIT:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
end

%% Inputs
if ~ischar(files_in.In1{1})
    error('files in should be a char');
end

[path_nii,name_nii,ext_nii] = fileparts(char(files_in.In1{1}));
if ~strcmp(ext_nii, '.nii')
    error('First file need to be a .nii');
end


if isfield(opt,'threshold') && (~isnumeric(opt.threshold))
    opt.threshold = str2double(opt.threshold);
    if isnan(opt.threshold)
        disp('The threshold used was not a number')
        return
    end
end


%% Building default output names
if strcmp(opt.folder_out,'') % if the output folder is left empty, use the same folder as the input
    opt.folder_out = path_nii;
end

if isempty(files_out)
    files_out.In1 = {cat(2,opt.folder_out,filesep,name_nii,'_',opt.T1_output_filename,ext_nii)};
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
N = niftiread(files_in.In1{1});
% load nifti info
info = niftiinfo(files_in.In1{1});

%% load input JSON file
J = spm_jsonread(strrep(files_in.In1{1}, '.nii', '.json'));

% retieve calculation parameter
first_scan=opt.first_scan;
M0mapYN=opt.M0_map;
IEmapYN=opt.IE_map;


%% retreive FAIR method
%T1map_method = scan_acqp('##$Method=',data.texte,0);
T1map_method = J.SequenceName.value{1};
if(~isempty(regexpi(T1map_method,'\w*fair\w*')))
    %PvVersion = scan_acqp('##TITLE=Parameter List, ParaVision ',data.texte,0);
    PvVersion = J.SoftwareVersions.value{1};
    if ~isnan(PvVersion(1)) && ~isempty(regexpi(PvVersion,'\w*6.\w*'))
        %FairMode=scan_acqp('##$PVM_FairMode=',data.texte,0); % SELECTIVE, NONSELECTIVE, INTERLEAVED or INTERLEAVED2
        FairMode = J.FairModePVM.value{1}; % SELECTIVE, NONSELECTIVE, INTERLEAVED or INTERLEAVED2
    else
        %FairMode=scan_acqp('##$FairMode=',data.texte,0); % SELECTIVE, NONSELECTIVE, INTERLEAVED or INTERLEAVED2
        FairMode = J.FairMode.value{1};
    end
    if(max(strcmp(FairMode,{'SELECTIVE','NONSELECTIVE'}))==1)
        NumFairMode = 1;
    elseif(max(strcmp(FairMode,{'INTERLEAVED','INTERLEAVED2'}))==1)
        NumFairMode = 2;
    end
else
    error('Wrong scan selected ?')
end

%% slice acquisition order for correct inversion time
%ObjOrderListName = sprintf('##$PVM_ObjOrderList=( %d )',size(data.reco.data,4));
if isfield(J, 'ObjOrderList')
    PVM_ObjOrderList = J.ObjOrderList.value;
end
%PVM_ObjOrderList=scan_acqp(ObjOrderListName,data.texte,1)';
if(size(N,3)==1)
    PVM_ObjOrderList = 0;
end

%% retrieve inversion time and inter slice time for different acquisition
% module of FAIR method (EPI, RARE, fisp, segm...)
if ( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ~isempty(regexpi(T1map_method,'\w*epi\w*')) )
    % inter slice time
    interSliceTime = InterSliceTimeFairEpi(J);
    
    % inversion time
    %FairTIR_Arr=(scan_acqp('FairTIR_Arr=',data.texte,1))';
    if ~isfield(J, 'FairTIRArr')
        FairTIR_Arr = J.FairTIRArrPVM.value.';
    else
        FairTIR_Arr = J.FairTIRArr.value.';
    end
    InvTimeRaw = FairTIR_Arr(first_scan:end);   % liste  des temps inversion pour la methode fair
    InvTimeRaw = repmat(InvTimeRaw,NumFairMode,1);            % premiere ligne selective deuxieme nonselective si les deux mesures sont faites
elseif( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ( ~isempty(regexpi(T1map_method,'\w*segm\w*')) ||...
        ~isempty(regexpi(T1map_method,'\w*fisp\w*')) ) )
    % inversion time
    %InvTimeSel = (scan_acqp('##$InvTimeSel=',data.texte,1))';
    InvTimeSel = J.InvTimeSel.value;
    %InvTimeGlo=(scan_acqp('##$InvTimeGlob=',data.texte,1))';
    InvTimeGlo = J.InvTimeGlo.value;
    if(strcmp(FairMode,'SELECTIVE'))
        InvTimeRaw=InvTimeSel(first_scan:end);
    elseif(strcmp(FairMode,'NONSELECTIVE'))
        InvTimeRaw=InvTimeGlo(first_scan:end);
    elseif(strcmp(FairMode,'INTERLEAVED'))
        InvTimeRaw(1,:)=InvTimeSel(first_scan:end);
        InvTimeRaw(2,:)=InvTimeGlo(first_scan:end);
    end
    % inter slice time
    interSliceTime = 0; % ms
    fprintf('Interlice time correction by default equal to %ims\n',interSliceTime)
    
    % T1map calculation based on a non selective FAIR
elseif ( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ~isempty(regexpi(T1map_method,'\w*rare\w*')) )
    % inversion time
    if isfield(J,'FairTIRArr')
        FairTIR_Arr = J.FairTIRArr.value.';
    elseif isfield(J,'FairTIRArrPVM')
        FairTIR_Arr = J.FairTIRArrPVM.value.';
    end
    InvTimeRaw = FairTIR_Arr(first_scan:end);   % liste  des temps inversion pour la methode fair
    InvTimeRaw = repmat(InvTimeRaw,NumFairMode,1);            % premiere ligne selective deuxieme nonselective si les deux mesures sont faites
    % inter slice time IS THIS CORRECT FOR RARE???
    interSliceTime = 0; % ms
    fprintf('Interlice time correction by default equal to %ims\n',interSliceTime)
end


% to allow single slice datasets to be processed
if size(size(N),2) > 3
    numberOfSlices = size(N,3); % To be checked on 4D or 5D datasets
else %single slice dataset
    numberOfSlices = 1;
end


if(NumFairMode==1)
    data_in_vector = reshape(N, [size(N,1)*size(N,2),...
        size(N,3), size(N,4)]);
elseif(NumFairMode==2)
    data_in_vector = permute(reshape(N, [size(N,1)*size(N,2),...
        size(N,3),NumFairMode, size(N,5)]), [1 2 4 3]);
end

% 2 cases bases on single slice or multislices
if numberOfSlices >1
    fit_T1_result = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
    if strcmp(M0mapYN,'Yes')
        fit_M0_result = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
        fit_M0_err = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
    end
    if strcmp(IEmapYN,'Yes')
        fit_IE_result = NaN(size(data_in_vector,1),size(data_in_vector,2),NumFairMode);
    end
else
    fit_T1_result = NaN(size(data_in_vector,1),NumFairMode);
    if strcmp(M0mapYN,'Yes')
        fit_M0_result = NaN(size(data_in_vector,1),NumFairMode);
    end
    if strcmp(IEmapYN,'Yes')
        fit_IE_result = NaN(size(data_in_vector,1),NumFairMode);
    end
end

%% Fitting

for it_mode = 1 : NumFairMode
    for it_slice = 1 : numberOfSlices
        data_in_vectorTemp = squeeze(data_in_vector(:,it_slice,:,it_mode));
        NumSlice = find(PVM_ObjOrderList==it_slice-1)-1;
        
        % retrieve parameters
        echotime_used = InvTimeRaw(it_mode,:) + NumSlice * interSliceTime;
        
        for voxel_nbr=1:size(data_in_vectorTemp,1)
            if sum(isnan(data_in_vectorTemp(voxel_nbr,:))) ~= length(data_in_vectorTemp(voxel_nbr,:))
                [~,voxel_min_nbr]=min(data_in_vectorTemp(voxel_nbr,:));
                [aaa, ~,  ccc]=levenbergmarquardt('fit_T1_3param',echotime_used,data_in_vectorTemp(voxel_nbr,:),[echotime_used(voxel_min_nbr)/0.693*1.2 max(data_in_vectorTemp(voxel_nbr,:)) 1]);
                if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 & ccc==-1 %#ok<AND2>
                    fit_T1_result(voxel_nbr,it_slice,it_mode)=aaa(1);
                    if strcmp(M0mapYN,'Yes')
                        fit_M0_result(voxel_nbr,it_slice,it_mode)=aaa(2);
                    end
                    if strcmp(IEmapYN,'Yes')
                        fit_IE_result(voxel_nbr,it_slice,it_mode)=aaa(3);
                    end
                end
            end
        end
        
    end
end

%% reshape and formatting of data
T1map=reshape(fit_T1_result,[size(N,1),size(N,2),numberOfSlices ,NumFairMode]); %

if strcmp(M0mapYN,'Yes')
    M0map=reshape(fit_M0_result,[size(N,1),size(N,2),numberOfSlices ,NumFairMode]); % LH20210120
end
if strcmp(IEmapYN,'Yes')
    IEmap=reshape(fit_IE_result,[size(N,1),size(N,2),numberOfSlices ,NumFairMode]); % LH20210120
end

%T1map.reco.data=permute(tmp, [1 2 4 3]);
% correction Kober 2004 on apparent T1
if ( ~isempty(regexpi(T1map_method,'\w*fair\w*')) &&...
        ( ~isempty(regexpi(T1map_method,'\w*segm\w*')) ||...
        ~isempty(regexpi(T1map_method,'\w*fisp\w*')) ) )
    PVM_ExcPulseAngle=scan_acqp('##$PVM_ExcPulseAngle=',data.texte,1);
    Seg_time=scan_acqp('##$Seg_time=',data.texte,1);
    %    T1map.reco.data = T1map.reco.data ./ ( 1 + T1map.reco.data * log( cos( PVM_ExcPulseAngle * pi / 180 ) ) / Seg_time );
end


OutputImages = T1map;

% save the new files (.nii & .json)
% update the header before saving the new .nii
info2 = info;
info2.Filename = files_out.In1{1};
info2.Filemoddate = char(datetime('now'));
info2.Datatype = class(OutputImages);
info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
info2.ImageSize = size(OutputImages);
info2.Description = [info.Description, 'Modified by Module_T1map_MIT Module'];

OutputImages(OutputImages < 0) = NaN;
OutputImages(OutputImages > 8000) = NaN;


% save the new .nii file
niftiwrite(OutputImages, files_out.In1{1}, info2);


%% Json Processing
[path, name, ~] = fileparts(files_in.In1{1});
jsonfile = [path, '/', name, '.json'];
J = ReadJson(jsonfile);
J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);

[path, name, ~] = fileparts(files_out.In1{1});
jsonfile = [path, '/', name, '.json'];
WriteJson(J, jsonfile)

% ESSAI PAS CONCLUANT - LH20210120 added for M0 and IE output
if strcmp(M0mapYN,'Yes')
    OutputImages = M0map;
    % save the new files (.nii & .json)
    % update the header before saving the new .nii
    info2 = info;
    info2.Filename = files_out.In2{1};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(OutputImages);
    info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
    info2.ImageSize = size(OutputImages);
    info2.Description = [info.Description, 'Modified by Module_T1map_MIT Module'];
    
    % save the new .nii file
    niftiwrite(OutputImages, files_out.In2{1}, info2);
    
    
    %% Json Processing
    [path, name, ~] = fileparts(files_in.In1{1});
    jsonfile = [path, '/', name, '.json'];
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    
    [path, name, ~] = fileparts(files_out.In2{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
    
end
if strcmp(IEmapYN,'Yes')
    OutputImages = IEmap;
    % save the new files (.nii & .json)
    % update the header before saving the new .nii
    info2 = info;
    info2.Filename = files_out.In3{1};
    info2.Filemoddate = char(datetime('now'));
    info2.Datatype = class(OutputImages);
    info2.PixelDimensions = info.PixelDimensions(1:length(size(OutputImages)));
    info2.ImageSize = size(OutputImages);
    info2.Description = [info.Description, 'Modified by Module_T1map_MIT Module'];
    
    % save the new .nii file
    niftiwrite(OutputImages, files_out.In3{1}, info2);
    
    
    %% Json Processing
    [path, name, ~] = fileparts(files_in.In1{1});
    jsonfile = [path, '/', name, '.json'];
    J = ReadJson(jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    
    [path, name, ~] = fileparts(files_out.In3{1});
    jsonfile = [path, '/', name, '.json'];
    WriteJson(J, jsonfile)
end
% LH20210120 added for M0 and IE output


%
end