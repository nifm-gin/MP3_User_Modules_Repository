function [files_in,files_out,opt] = Module_SPM_Dartel_create_Templates(files_in,files_out,opt)

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
    module_option(:,3)   = {'Flowfield_filename_prefix','u_'};
    module_option(:,4)   = {'output_template_Name','Template'};
    module_option(:,5)   = {'Regularisation_Form','Linear Elastic Energy'};
    module_option(:,6)   = {'Outer_iterations', 6};
    module_option(:,7)   = {'Inner_Iterations', '3; 3; 3; 3; 3; 3'};
    module_option(:,8)   = {'Reg_params', '4 2 1e-06; 2 1 1e-06; 1 0.5 1e-06; 0.5 0.25 1e-06; 0.25 0.125 1e-06; 0.25 0.125 1e-06'};
    module_option(:,9)   = {'Time_Steps', '1; 1; 2; 4; 16; 64'};
    module_option(:,10)   = {'Smoothing_Parameter', '16; 8; 4; 2; 1; 0.5'};
    module_option(:,11)   = {'LM_Regularisation', 0.01};
    module_option(:,12)   = {'Cycle', 3};
    module_option(:,13)   = {'Iterations', 3};
    module_option(:,14)   = {'AutomaticJobsCreation', 'No'};
    module_option(:,15)   = {'RefInput',1};
    module_option(:,16)   = {'InputToReshape',1};
    module_option(:,17)   = {'Table_in', table()};
    module_option(:,18)   = {'Table_out', table()};
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
    'Run the Dartel nonlinear image registration procedure.'
    'This involves iteratively matching all the selected images to a template generated from their own mean. A series of Template*.nii files are generated, which become increasingly crisp as the registration proceeds.'
    ''
    'This module will generate :'
    '   - serveral templates based on the number of iterations asked by the user. So far, these templates are not saved in the database (but will be generated in the /tmp folder'
    '   - Flow field maps for each inputs that parameterise the deformations for each input to the template generated during this module'
    }'};

    user_parameter(:,2)   = {'Input Scans','XScan','','', {'SequenceName'},'Mandatory',...
         {'Select the images to be warped together. Multiple sets of images can be simultaneously registered.'
         'For example, used can use the grey and white matter images.'
         }'};
    user_parameter(:,3)   = {'Parameters','','','','','','Various settings for the optimisation. The default values should work reasonably well for aligning tissue class images together.'};
    user_parameter(:,4)   = {'   .Template basename','char','','output_template_Name','','',...
        'This module will create serveral template based on the number of iterations asked by the user. So far, these templates are not saved in the database (but will be generated in the /tmp folder'};
    user_parameter(:,5)  = {'   .Flowfield filename prefix','char','','Flowfield_filename_prefix','','',...
        'Specify the string to be prepended to the filenames of the flowfield image. Default prefix is ''u_''.'};
    user_parameter(:,6)   = {'   .Regularisation Form','cell', {'Linear Elastic Energy', 'Membrane Energy', 'Bending Energy'},'Regularisation_Form','','',...
        'The registration is penalised by some ``energy'''' term.  Here, the form of this energy term is specified. Three different forms of regularisation can currently be used.'};
    user_parameter(:,7)   = {'   .Outer Iterations','','','','','','The images are averaged, and each individual image is warped to match this average. This is repeated a number of times.'};
    user_parameter(:,8)   = {'      .Number of outer iteration','numeric','','Outer_iterations','','',...
        {'Outer Iteration Different parameters can be specified for each outer iteration.'
        'Each of the warps the images to the template, and then regenerates the template from the average of the warped images.'
        'Multiple outer iterations should be used for more accurate results, beginning with a more coarse registration (more regularisation) then ending with the more detailed registration (less regularisation).'
        }'};
    
    user_parameter(:,9)   = {'      .Inner Iterations','numeric','','Inner_Iterations','','',...
        {'The number of Gauss-Newton iterations to be done within this outer iteration.'
        'After this, new average(s) are created, which the individual images are warped to match.'
        }' };
    user_parameter(:,10)   = {'      .Reg params','numeric','','Reg_params','','',...
        {'For linear elasticity, the parameters are mu, lambda and id.'
        'For membrane energy, the parameters are lambda, unused and id.id is a term for penalising absolute displacements, and should therefore be small.' 
        'For bending energy, the parameters are lambda, id1 and id2, and the regularisation is by (-lambda*Laplacian + id1)?2 + id2.'
        ''
        'Use more regularisation for the early iterations so that the deformations are smooth, and then use less for the later ones so that the details can be better matched.'
        ''
        'An Number_of_outer_iteration time 1-by-3 vector separated by a ";" must be entered'}' };
   user_parameter(:,11)   = {'      .Time Steps','numeric','','Time_Steps','','',...
        {'The number of time points used for solving the partial differential equations.'
        'A single time point would be equivalent to a small deformation model.'
        'Smaller values allow faster computations, but are less accurate in terms of inverse consistency and may result in the one-to-one mapping breaking down.'
        'Earlier iteration could use fewer time points, but later ones should use about 64 (or fewer if the deformations are very smooth).'
        ''
        'An Number_of_outer_iteration time a scalar separated by a ";" must be entered'}' };
    user_parameter(:,12)   = {'      .Smoothing Parameter','numeric','','Smoothing_Parameter','','',...
        {'A LogOdds parameterisation of the template is smoothed using a multi-grid scheme. The amount of smoothing is determined by this parameter.'
        ''
        'An Number_of_outer_iteration time a scalar separated by a ";" must be entered'}' };   
    user_parameter(:,13)   = {'   .Optimisation Settings','','','','','',...
        {'Settings for the optimisation. If you are unsure about them, then leave them at the default values.' 
        'Optimisation is by repeating a number of Levenberg-Marquardt iterations, in which the equations are solved using a full multi-grid (FMG) scheme.'
        'FMG and Levenberg-Marquardt are both described in Numerical Recipes (2nd edition).'
        }'};   
    user_parameter(:,14)   = {'      .LM Regularisation','numeric','','LM_Regularisation','','',...
        {'Levenberg-Marquardt regularisation. Larger values increase the the stability of the optimisation, but slow it down.'
        'A value of zero results in a Gauss-Newton strategy, but this is not recommended as it may result in instabilities in the FMG.'
        ''
        'An 1-by-1 array must be entered'}' };
    user_parameter(:,15)   = {'      .Cycle','numeric','','Cycle','','',...
        {'Number of cycles used by the full multi-grid matrix solver.'
        'More cycles result in higher accuracy, but slow down the algorithm.'
        'See Numerical Recipes for more information on multi-grid methods.'
        ''
        'Please enter a value between 1 to 8; default = 3'}' };   
    user_parameter(:,16)   = {'      .Iterations','numeric','','Iterations','','',...
        {'Number of relaxation iterations performed in each multi-grid cycle.'
        'More iterations are needed if using ?bending energy? regularisation, because the relaxation scheme only runs very slowly.'
        'See the chapter on solving partial differential equations in Numerical Recipes for more information about relaxation methods.'
        ''
        'Please enter a value between 1 to 8; default = 3'}' };   


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
% check if user enter correctly all input parameters needed
if numel(str2num(opt.Inner_Iterations)) ~= opt.Outer_iterations
    error('Module_SPM_Dartel_create_Templates:brick','Number of Inner_Iterations does not correspond to the number of outer interation, type ''help %s'' for more info.',mfilename)
elseif isempty(str2num(opt.Reg_params)) || size(str2num(opt.Reg_params),1) ~= opt.Outer_iterations
    error('Module_SPM_Dartel_create_Templates:brick','Number Reg parameters (or it format n time 1-by-3) does not correspond to the number of outer interation, type ''help %s'' for more info.',mfilename)
elseif numel(str2num(opt.Time_Steps)) ~= opt.Outer_iterations
    error('Module_SPM_Dartel_create_Templates:brick','Number of time steps does not correspond to the number of outer interation, type ''help %s'' for more info.',mfilename)
elseif numel(str2num(opt.Smoothing_Parameter)) ~= opt.Outer_iterations
    error('Module_SPM_Dartel_create_Templates:brick','Number of smoothing parameters does not correspond to the number of outer interation, type ''help %s'' for more info.',mfilename)
end





%%%
if strcmp(files_out, '')
    % create outputs for the templates
    % this will be coded later. First we need to find a way to figure out
    % where the template(s) will be saved (which patient/timepoint)
        
    
    % only patient/timepoint containing every inputs will be processed
    UniqueTable = table;
    UniqueTable.Patient = opt.Table_in.Patient;
    UniqueTable.Tp = opt.Table_in.Tp;
    UniqueTable = unique(UniqueTable, 'rows');
    % check every patient/timepoints
    for i = 1:size(UniqueTable,1)
        % if the current patient/timepoint contain all scans required -->
        % create an ouput scan/table
        if isequal(unique(sort(opt.Table_in.SequenceName((opt.Table_in.Patient == UniqueTable.Patient(i)) & (opt.Table_in.Tp == UniqueTable.Tp(i))))), sort(unique(opt.Table_in.SequenceName)))
            tmp_current_table_out = opt.Table_in(opt.Table_in.Patient == UniqueTable.Patient(i) & opt.Table_in.Tp == UniqueTable.Tp(i),:);
            % use the first scan as model
            tmp_current_table_out = tmp_current_table_out(1,:);
            tmp_current_table_out.IsRaw = categorical(0);
            tmp_current_table_out.Path = categorical(cellstr([opt.folder_out, filesep]));
            tmp_current_table_out.SequenceName = categorical(cellstr([opt.Flowfield_filename_prefix, char(tmp_current_table_out.SequenceName)]));
            tmp_current_table_out.Filename = categorical(cellstr([char(tmp_current_table_out.Patient), '_', char(tmp_current_table_out.Tp), '_', char(tmp_current_table_out.SequenceName)]));
            f_out = [char(tmp_current_table_out.Path), char(tmp_current_table_out.Patient), '_', char(tmp_current_table_out.Tp), '_', char(tmp_current_table_out.SequenceName), '.nii'];
            if  ~isfield(files_out, 'In2')
                files_out.In2{1} = f_out;
            else
                files_out.In2{length(files_out.In2)+1} = f_out;
            end
            opt.Table_out =  [opt.Table_out; tmp_current_table_out];
        end
    end
    
end




%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_SPM_Dartel_create_Templates:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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


input_scans_name = unique(opt.Table_in.SequenceName);
% for i = 1:numel(inputs_listing)
%     matlabbatch{1}.spm.tools.dartel.warp.images = {{i}}
% end



for  i = 1:numel(input_scans_name)
    inputs = cell(1,size(opt.Table_out,1));
    for j=1:size(opt.Table_out,1)
        % first, copy all inputs to the /tmp fold
        % Otherwise, spm will generate the output scans in the wrong
        % directory
        curent_Table = opt.Table_in(opt.Table_in.Patient == opt.Table_out.Patient(j) & opt.Table_in.Tp == opt.Table_out.Tp(j) & opt.Table_in.SequenceName == input_scans_name(i), :);
        curent_nii_filename = [char(curent_Table.Path) char(curent_Table.Filename) '.nii'];
        if ~strcmp(curent_nii_filename,  strrep(curent_nii_filename, '/Derived_data/', '/Tmp/'))
            copyfile(curent_nii_filename,  strrep(curent_nii_filename, '/Derived_data/', '/Tmp/'));
        end
        % save the inputs in a SPM format
        inputs{j} = [strrep(curent_nii_filename, '/Derived_data/', '/Tmp/'), ',1'];   
    end
    matlabbatch{1}.spm.tools.dartel.warp.images{i} = inputs';
end

matlabbatch{1}.spm.tools.dartel.warp.settings.template = opt.output_template_Name;
switch opt.Regularisation_Form
    case 'Linear Elastic Energy'
        matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
    case 'Membrane Energy'
        matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 1;
    case 'Bending Energy'
        matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 1;
end

its = str2num(opt.Inner_Iterations); %#ok<*ST2NM>
rparam = str2num(opt.Reg_params);
K = str2num(opt.Time_Steps);
slam = str2num(opt.Smoothing_Parameter);

for i = 1:opt.Outer_iterations
    matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).its = its(i);
    matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).rparam =  rparam(i,:);
    switch K(i)
        case 1
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 0;
        case 2
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 1;
        case 4
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 2;
        case 8
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 3;
        case 16
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 4;
        case 32
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 5;
        case 64
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 6;
        case    '128'
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 7;
        case    '256'
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 8;
        case    '512'
            matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).K = 9;
    end
    matlabbatch{1}.spm.tools.dartel.warp.settings.param(i).slam = slam(i);
    
end

matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = opt.LM_Regularisation;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = opt.Cycle;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = opt.Iterations;


jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
tic
spm_jobman('run', jobs, inputs{:});
toc


% rename output files in order to match to MP3 opt.Table_out
for j=1:size(opt.Table_out,1)
    % copy the json of the first input
    [path, name, ext] = fileparts(strrep(matlabbatch{1}.spm.tools.dartel.warp.images{1}{j}, ',1', ''));
    curent_Table = opt.Table_in(opt.Table_in.Filename == name,:);

    %curent_Table = opt.Table_in(opt.Table_in.Patient == opt.Table_out.Patient(j) & opt.Table_in.Tp == opt.Table_out.Tp(j) & opt.Table_in.SequenceName == input_scans_name(i), :);
    curent_json_filename = [char(curent_Table.Path) char(curent_Table.Filename) '.json'];
    new_jsonfile = [char(opt.Table_out.Path(j)) char(opt.Table_out.Filename(j)) '.json'];

    copyfile(curent_json_filename, new_jsonfile )
    %% Json Processing
%     [path, name, ~] = fileparts(files_in.In1{i});
    J = ReadJson(new_jsonfile);
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    WriteJson(J, new_jsonfile)
    
    % rename the nii file
    %spm_filename = strrep(matlabbatch{1}.spm.tools.dartel.warp.images{1}{j}, '.nii,1', '.nii');
   % [path, name, ext] = fileparts(spm_filename);
    movefile(fullfile(path, [opt.Flowfield_filename_prefix name '_' opt.output_template_Name ext]), files_out.In2{j})
       
end
    
    





