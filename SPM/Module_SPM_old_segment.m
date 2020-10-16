function [files_in,files_out,opt] = Module_SPM_old_segment(files_in,files_out,opt)

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
    module_option(:,3)   = {'Ouput_grey_Matter_opt','Native Space'};
    module_option(:,4)   = {'Ouput_white_Matter_opt','Native Space'};
    module_option(:,5)   = {'Ouput_csf_Matter_opt','Native Space'};
    module_option(:,6)   = {'Bias_corrected','Save Bias Corrected'};
    module_option(:,7)   = {'Clean_up_partition','Dont do cleanup'};
    module_option(:,8)   = {'Gaussians_per_class','2 2 2 4'};
    module_option(:,9)   = {'Affine_regularisation','ICBM space template - European brains'};
    module_option(:,10)   = {'Warping_regularisation','1'};
    module_option(:,11)   = {'Warp_frequency_cutoff','25'};
    module_option(:,12)   = {'Bias_regularisation','very light regularisation (0.0001)'};
    module_option(:,13)   = {'Bias_FWHM','60mm cutoff'};
    module_option(:,14)   = {'Sampling_distance','3'};

    
    module_option(:,15)   = {'RefInput',1};
    module_option(:,16)   = {'InputToReshape',[2 3 4 5]};
    
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
        'This uses a crude routine for extracting the brain from segmented images.  It begins by taking the white matter, and eroding it acouple of times to get rid of any odd voxels.  The algorithm continues on to do conditional dilations for several iterations, where the condition is based upon gray or white matter being present. This identified region is then used to clean up the grey and white matter partitions, and has a slight influences on the CSF partition.'
        ''
        'If you find pieces of brain being chopped out in your data, then you may wish to disable or tone down the cleanup procedure.'
        ''
        'This routine produces spatial normalisation parameters (*_seg_sn.mat files) by default. These can be used for writing spatially normalised versions of your data, via the "Normalise: Write" option. This mechanism may produce superior results than the "Normalise: Estimate" option (but probably not as good as those produced using DARTEL).'
        ''
        'In addition, it also produces files that can be used for doing inverse normalisation. If you have an image of regions defined in the standard space, then the inverse deformations can be used to warp these regions so that it approximately overlay your image. To use this facility, the bounding-box and voxel sizes should be set to non-finite values (e.g. [NaN NaN NaN] for the voxel sizes, and ones(2,3)*NaN for the bounding box. This would be done by the spatial normalisation module, which allows you to select a set of parameters that describe the nonlinear warps, and the images that they should be applied to.'
        ''
        'There are a number of options about what data you would like the routine to produce. The routine can be used for producing images of tissue classes, as well as bias corrected images. The native space option will produce a tissue class image (c*) that is in alignment with the original/* (see Figure \ref{seg1})*/.  You can also produce spatially normalised versions - both with (mwc*) and without (wc*) modulation/* (see Figure \ref{seg2})*/. The bounding box and voxel sizes of the spatially normalised versions are the same as that of the tissue probability maps with which they are registered. These can be used for doing voxel-based morphometry with (also see the ``Using DARTEL'' chapter of the manual). All you need to do is smooth them and do the stats (which means no more questions on the mailing list about how to do "optimized VBM").'
        ''
        'Modulation is to compensate for the effect of spatial normalisation.  When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images.  For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labelled grey matter.  In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping.  If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.'
        '/*\begin{figure} \begin{center} \includegraphics[width=140mm]{images/seg1} \end{center} \caption{Segmentation results. These are the results that can be obtained in the original space of the image (i.e. the results that are not spatially normalised). Top left: original image (X.img). Top right: bias corrected image (mX.img). Middle and bottom rows: segmented grey matter (c1X.img), white matter (c2X.img) and CSF (c3X.img). \label{seg1}} \end{figure} */'
        '/*\begin{figure} \begin{center} \includegraphics[width=140mm]{images/seg2} \end{center} \caption{Segmentation results. These are the spatially normalised results that can be obtained (note that CSF data is not shown). Top row: The tissue probability maps used to guide the segmentation. Middle row: Spatially normalised tissue maps of grey and white matter (wc1X.img and wc2X.img). Bottom row: Modulated spatially normalised tissue maps of grey and white matter (mwc1X.img and mwc2X.img). \label{seg2}} \end{figure} */'
        'A deformation field is a vector field, where three values are associated with each location in the field.  The field maps from co-ordinates in the normalised image back to co-ordinates in the original image.  The value of the field at co-ordinate [x y z] in the normalised space will be the co-ordinate [x'' y'' z''] in the original volume. The gradient of the deformation field at a co-ordinate is its Jacobian matrix, and it consists of a 3x3 matrix:'
        ''
        '%   /                      \'
        '%   | dx''/dx  dx''/dy dx''/dz |'
        '%   |                       |'
        '%   | dy''/dx  dy''/dy dy''/dz |'
        '%   |                       |'
        '%   | dz''/dx  dz''/dy dz''/dz |'
        '%   \                      /'
        '/* \begin{eqnarray*}\begin{pmatrix}\frac{dx''}{dx} & \frac{dx''}{dy} & \frac{dx''}{dz}\cr\frac{dy''}{dx} & \frac{dy''}{dy} & \frac{dy''}{dz}\cr\frac{dz''}{dx} & \frac{dz''}{dy} & \frac{dz''}{dz}\cr\end{pmatrix}\end{eqnarray*}*/'
        'The value of dx''/dy is a measure of how much x'' changes if y is changed by a tiny amount. The determinant of the Jacobian is the measure of relative volumes of warped and unwarped structures.'
        'The modulation step simply involves multiplying by the relative volumes /*(see Figure \ref{seg2})*/.'
        }'};
    user_parameter(:,2)   = {'Reference Image','1Scan','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
        'Select scans for processing. This assumes that there is one scan for each subject. Note that multi-spectral (when there are two or more registered images of different contrasts) processing is not yet implemented for this method.'};
    user_parameter(:,3)   = {'Output Files','','','','','',''};
    user_parameter(:,4)   = {'    .Grey Matter','cell', {'None','Native Space', 'Unmodulated Normalised', 'Modulated Normalised', 'Native + Unmodulated Normalised', 'Native + Modulated Normalised','Native + Modulated + Unmodulated', 'Modulated + Unmodulated Normalised'},'Ouput_grey_Matter_opt','','',...
        {
        'Options to produce grey matter images: c1*, wc1* and mwc1*.'
        'None'
        'Native Space'
        'Unmodulated Normalised'
        'Modulated Normalised'
        'Native + Unmodulated Normalised'
        'Native + Modulated Normalised'
        'Native + Modulated + Unmodulated'
        'Modulated + Unmodulated Normalised'
        }'};
    user_parameter(:,5)   = {'    .White Matter','cell', {'None','Native Space', 'Unmodulated Normalised', 'Modulated Normalised', 'Native + Unmodulated Normalised', 'Native + Modulated Normalised','Native + Modulated + Unmodulated', 'Modulated + Unmodulated Normalised'},'Ouput_white_Matter_opt','','',...
        {
        'Options to produce White matter images: c1*, wc1* and mwc1*.'
        'None'
        'Native Space'
        'Unmodulated Normalised'
        'Modulated Normalised'
        'Native + Unmodulated Normalised'
        'Native + Modulated Normalised'
        'Native + Modulated + Unmodulated'
        'Modulated + Unmodulated Normalised'
        }'};
    
    user_parameter(:,6)   = {'    .Cerebro-Spinal Fluid','cell', {'None','Native Space', 'Unmodulated Normalised', 'Modulated Normalised', 'Native + Unmodulated Normalised', 'Native + Modulated Normalised','Native + Modulated + Unmodulated', 'Modulated + Unmodulated Normalised'},'Ouput_csf_Matter_opt','','',...
        {
        'Options to produce Cerebro-Spinal Fluid images: c1*, wc1* and mwc1*.'
        'None'
        'Native Space'
        'Unmodulated Normalised'
        'Modulated Normalised'
        'Native + Unmodulated Normalised'
        'Native + Modulated Normalised'
        'Native + Modulated + Unmodulated'
        'Modulated + Unmodulated Normalised'
        }'};
    user_parameter(:,7)  =  {'    .Bias corrected','cell', {'Save Bias Corrected', 'Don''t Save Corrected'},'Bias_corrected','','',...
        'This is the option to produce a bias corrected version of your image. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'};
    user_parameter(:,8)  =  {'    .Clean up and partition','cell',{'Dont do cleanup', 'Light Clean', 'Thorough Clean'},'Clean_up_partition','','',...
        {'This uses a crude routine for extracting the brain from segmented images.  It begins by taking the white matter, and eroding it acouple of times to get rid of any odd voxels.  The algorithm continues on to do conditional dilations for several iterations, where the condition is based upon gray or white matter being present. This identified region is then used to clean up the grey and white matter partitions, and has a slight influences on the CSF partition.'
        ''
        'If you find pieces of brain being chopped out in your data, then you may wish to disable or tone down the cleanup procedure.'}'
        };
    user_parameter(:,9)  = {'Custom','','','','','',...
        'Various options can be adjusted in order to improve the performance of the algorithm with your data.  Knowing what works best should be a matter of empirical exploration.  For example, if your data has very little intensity non-uniformity artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'};
    user_parameter(:,10)  = {'    .Tissue probability maps','','','','','',...
        {'Select the tissue probability images. These should be maps of grey matter, white matter and cerebro-spinal fluid probability. A nonlinear deformation field is estimated that best overlays the tissue probability maps on the individual subjects'' image. The default tissue probability maps are modified versions of the ICBM Tissue Probabilistic Atlases.These tissue probability maps are kindly provided by the International Consortium for Brain Mapping, John C. Mazziotta and Arthur W. Toga. http://www.loni.ucla.edu/ICBM/ICBM_TissueProb.html. The original data are derived from 452 T1-weighted scans, which were aligned with an atlas space, corrected for scan inhomogeneities, and classified into grey matter, white matter and cerebrospinal fluid. These data were then affine registered to the MNI space and downsampled to 2mm resolution.'
        ''
        'Rather than assuming stationary prior probabilities based upon mixing proportions, additional information is used, based on other subjects'' brain images.  Priors are usually generated by registering a large number of subjects together, assigning voxels to different tissue types and averaging tissue classes over subjects. Three tissue classes are used: grey matter, white matter and cerebro-spinal fluid. A fourth class is also used, which is simply one minus the sum of the first three. These maps give the prior probability of any voxel in a registered image being of any of the tissue classes - irrespective of its intensity.'
        ''
        'The model is refined further by allowing the tissue probability maps to be deformed according to a set of estimated parameters. This allows spatial normalisation and segmentation to be combined into the same model. This implementation uses a low-dimensional approach, which parameterises the deformations by a linear combination of about a thousand cosine transform bases. This is not an especially precise way of encoding deformations, but it can model the variability of overall brain shape. Evaluations by Hellier et al have shown that this simple model can achieve a registration accuracy comparable to other fully automated methods with many more parameters.'
        }'};
    user_parameter(:,11)  = {'        .Grey Matter probabiliy map','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
        'Please select the Grey Matter probability map'};
    user_parameter(:,12)  = {'        .White Matter probabiliy map','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
        'Please select the Grey Matter probability map'};
    user_parameter(:,13)  = {'        .Cerebro-Spinal Fluid probabiliy map','1Scan1TPXP','','', {'SequenceName', 'Tp', 'Patient'},'Mandatory',...
        'Please select the Grey Matter probability map'};
    user_parameter(:,14)  = {'    .Gaussians per class','numeric','','Gaussians_per_class','','',...
        'The number of Gaussians used to represent the intensity distribution for each tissue class can be greater than one. In other words, a tissue probability map may be shared by several clusters. The assumption of a single Gaussian distribution for each class does not hold for a number of reasons. In particular, a voxel may not be purely of one tissue type, and instead contain signal from a number of different tissues (partial volume effects). Some partial volume voxels could fall at the interface between different classes, or they may fall in the middle of structures such as the thalamus, which may be considered as being either grey or white matter. Various other image segmentation approaches use additional clusters to model such partial volume effects. These generally assume that a pure tissue class has a Gaussian intensity distribution, whereas intensity distributions for partial volume voxels are broader, falling between the intensities of the pure classes. Unlike these partial volume segmentation approaches, the model adopted here simply assumes that the intensity distribution of each class may not be Gaussian, and assigns belonging probabilities according to these non-Gaussian distributions. Typical numbers of Gaussians could be two for grey matter, two for white matter, two for CSF, and four for everything else.'};
    user_parameter(:,15)  = {'    .Affine Regularisation','cell',{'No Affine Registration', 'ICBM space template - European brains', 'ICBM space template - East Asian brains', 'Average sized template', 'No regularisation'},'Affine_regularisation','','',...
        {'The procedure is a local optimisation, so it needs reasonable initial starting estimates. Images should be placed in approximate alignment using the Display function of SPM before beginning. A Mutual Information affine registration with the tissue probability maps (D''Agostino et al, 2004) is used to achieve approximate alignment. Note that this step does not include any model for intensity non-uniformity. This means that if the procedure is to be initialised with the affine registration, then the data should not be too corrupted with this artifact.If there is a lot of intensity non-uniformity, then manually position your image in order to achieve closer starting estimates, and turn off the affine registration.'
         ''
         'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). For example, if registering to an image in ICBM/MNI space, then choose this option.  If registering to a template that is close in size, then select the appropriate option for this.'
         }'};
    user_parameter(:,16)  = {'    .Warping Regularisation','numeric','','Warping_regularisation','','',...
        {'The objective function for registering the tissue probability maps to the image to process, involves minimising the sum of two terms. One term gives a function of how probable the data is given the warping parameters. The other is a function of how probable the parameters are, and provides a penalty for unlikely deformations. Smoother deformations are deemed to be more probable. The amount of regularisation determines the tradeoff between the terms. Pick a value around one.  However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude). More regularisation gives smoother deformations, where the smoothness measure is determined by the bending energy of the deformations. '
        }'};
    user_parameter(:,17)  = {'    .Warp Frequency Cutoff','numeric','','Warp_frequency_cutoff','','',...
        {'Cutoff of DCT bases.  Only DCT bases of periods longer than the cutoff are used to describe the warps. The number actually used will depend on the cutoff and the field of view of your image. A smaller cutoff frequency will allow more detailed deformations to be modelled, but unfortunately comes at a cost of greatly increasing the amount of memory needed, and the time taken.'
        }'};
    user_parameter(:,18)  = {'    .Bias regularisation','cell',{'no regularisation (0)', 'extremely light regularisation (0.00001)', 'very light regularisation (0.0001)', 'light regularisation (0.001)', 'medium regularisation (0.01)', 'heavy regularisation (0.1)', 'very heavy regularisation (1)', 'extremely heavy regularisation (10)'},'Bias_regularisation','','',...
        {'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
          ''
         'An important issue relates to the distinction between intensity variations that arise because of bias artifact due to the physics of MR scanning, and those that arise due to different tissue properties.  The objective is to model the latter by different tissue classes, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity non-uniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of non-uniformity.'
        }'};
    user_parameter(:,19)  = {'    .Bias FWHM','cell',{'30mm cutoff', '40mm cutoff', '50mm cutoff', '60mm cutoff', '70mm cutoff', '80mm cutoff', '90mm cutoff', '100mm cutoff', '110mm cutoff', '120mm cutoff', '130mm cutoff', '140mm cutoff', '150mm cutoff', 'No correction'},'Bias_FWHM','','',...
        {'FWHM of Gaussian smoothness of bias. If your intensity non-uniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity non-uniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities.'
        }'};
    user_parameter(:,20)  = {'    .Sampling distance','numeric','','Sampling_distance','','',...
        {'The approximate distance between sampled points when estimating the model parameters. Smaller values use more of the data, but the procedure is slower.'
        }'};
    user_parameter(:,21)  = {'    .Masking image','1ROI1TPXP','','',{'SequenceName', 'Tp', 'Patient'},'Optional',...
        {'The segmentation can be masked by an image that conforms to the same space as the images to be segmented.  If an image is selected, then it must match the image(s) voxel-for voxel, and have the same voxel-to-world mapping.  Regions containing a value of zero in this image do not contribute when estimating the various parameters. '
        }'};
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


if strcmp(files_out, '')
    [Path_In1, Name_In1, ~] = fileparts(files_in.In1{1});
    tags1 = opt.Table_in(opt.Table_in.Path == [Path_In1, filesep],:);
    tags1 = tags1(tags1.Filename == Name_In1,:);
    assert(size(tags1, 1) == 1);
    if isfield(files_in, 'In2') && ~isempty(files_in.In2) && ~strcmp(opt.Ouput_grey_Matter_opt, 'None')
        tags_out_tmp = tags1;
        tags_out_tmp.IsRaw = categorical(0);
        tags_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        tags_out_tmp.SequenceName = categorical(cellstr(['c1', char(tags_out_tmp.SequenceName)]));
        tags_out_tmp.Filename = categorical(cellstr([char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName)]));
        f_out = [char(tags_out_tmp.Path), char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName), '.nii'];
        files_out.In1{1} = f_out;
        opt.Table_out = [opt.Table_out ; tags_out_tmp];
    end
    if isfield(files_in, 'In3') && ~isempty(files_in.In3) && ~strcmp(opt.Ouput_white_Matter_opt, 'None')
        tags_out_tmp = tags1;
        tags_out_tmp.IsRaw = categorical(0);
        tags_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        tags_out_tmp.SequenceName = categorical(cellstr(['c2', char(tags_out_tmp.SequenceName)]));
        tags_out_tmp.Filename = categorical(cellstr([char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName)]));
        f_out = [char(tags_out_tmp.Path), char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName), '.nii'];
        files_out.In2{1} = f_out;
        opt.Table_out = [opt.Table_out ; tags_out_tmp];  
    end
    if isfield(files_in, 'In4') && ~isempty(files_in.In4) && ~strcmp(opt.Ouput_csf_Matter_opt, 'None')
        tags_out_tmp = tags1;
        tags_out_tmp.IsRaw = categorical(0);
        tags_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
        tags_out_tmp.SequenceName = categorical(cellstr(['c3', char(tags_out_tmp.SequenceName)]));
        tags_out_tmp.Filename = categorical(cellstr([char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName)]));
        f_out = [char(tags_out_tmp.Path), char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName), '.nii'];
        files_out.In3{1} = f_out;
        opt.Table_out = [opt.Table_out ; tags_out_tmp];     
    end
    % add the _seg_inv_sn.mat as output
    tags_out_tmp = tags1;
    tags_out_tmp.IsRaw = categorical(0);
    tags_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_tmp.Type = categorical(cellstr('Mfile'));
    tags_out_tmp.SequenceName = categorical(cellstr([char(tags_out_tmp.SequenceName), '_seg_inv_sn']));
    tags_out_tmp.Filename = categorical(cellstr([char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName)]));
    f_out = [char(tags_out_tmp.Path), char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName), '.mat'];
    files_out.In4{1} = f_out;
    opt.Table_out = [opt.Table_out ; tags_out_tmp];
    % add the _seg_sn.mat as output
    tags_out_tmp = tags1;
    tags_out_tmp.IsRaw = categorical(0);
    tags_out_tmp.Path = categorical(cellstr([opt.folder_out, filesep]));
    tags_out_tmp.Type = categorical(cellstr('Mfile'));
    tags_out_tmp.SequenceName = categorical(cellstr([char(tags_out_tmp.SequenceName), '_seg_sn']));
    tags_out_tmp.Filename = categorical(cellstr([char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName)]));
    f_out = [char(tags_out_tmp.Path), char(tags_out_tmp.Patient), '_', char(tags_out_tmp.Tp), '_', char(tags_out_tmp.SequenceName), '.mat'];
    files_out.In5{1} = f_out;
    opt.Table_out = [opt.Table_out ; tags_out_tmp];

end    





%% Syntax
if ~exist('files_in','var')||~exist('files_out','var')||~exist('opt','var')
    error('Module_SPM_old_segement:brick','Bad syntax, type ''help %s'' for more info.',mfilename)
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
%%

% % First duplicate the reference image to the tmp folder and update the
% files_in.In1 variable

new_pathname = strcat(char(opt.Table_out.Path(1)),   char(opt.Table_in.Filename(1)), '.nii');
% if this module is executed in a pipeline, is it possible that
% files_in.In1{1} is already in the Tmp --> if so, we rename it
if strcmp(files_in.In1{1}, new_pathname)
    new_pathname = strrep(new_pathname, char(opt.Table_in.Filename(1)), [char(opt.Table_in.Filename(1)), '-bis']);
    opt.Table_in.Filename(1) = categorical(cellstr([char(opt.Table_in.Filename(1)), '-bis']));
    opt.Table_in.SequenceName(1) = categorical(cellstr([char(opt.Table_in.SequenceName(1)), '-bis']));
end

copyfile(files_in.In1{1},   new_pathname)
copyfile(strrep(files_in.In1{1},'.nii','.json'),  strrep(new_pathname,'.nii','.json'))
files_in.In1{1} = new_pathname;


matlabbatch{1}.spm.tools.oldseg.data(1) = {[files_in.In1{1}, ',1']};
switch opt.Ouput_grey_Matter_opt
    case 'None'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 0];
    case 'Native Space'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 1];
    case 'Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 1 0];
    case 'Modulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [1 0 0];
    case 'Native + Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [0 1 1];
    case 'Native + Modulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [1 0 1];
    case 'Native + Modulated + Unmodulated'
        matlabbatch{1}.spm.tools.oldseg.output.GM = [1 1 1];
    case 'Modulated + Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.GM =  [1 1 0];
end
switch opt.Ouput_white_Matter_opt
    case 'None'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 0];
    case 'Native Space'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 1];
    case 'Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 1 0];
    case 'Modulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [1 0 0];
    case 'Native + Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [0 1 1];
    case 'Native + Modulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [1 0 1];
    case 'Native + Modulated + Unmodulated'
        matlabbatch{1}.spm.tools.oldseg.output.WM = [1 1 1];
    case 'Modulated + Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.WM =  [1 1 0];
end
switch opt.Ouput_csf_Matter_opt
    case 'None'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 0];
    case 'Native Space'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 1];
    case 'Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 1 0];
    case 'Modulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [1 0 0];
    case 'Native + Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 1 1];
    case 'Native + Modulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [1 0 1];
    case 'Native + Modulated + Unmodulated'
        matlabbatch{1}.spm.tools.oldseg.output.CSF = [1 1 1];
    case 'Modulated + Unmodulated Normalised'
        matlabbatch{1}.spm.tools.oldseg.output.CSF =  [1 1 0];
end
if strcmp(opt.Bias_corrected, 'Don''t Save Corrected')
    matlabbatch{1}.spm.tools.oldseg.output.biascor = 0;
else
    matlabbatch{1}.spm.tools.oldseg.output.biascor = 1;
end
if strcmp(opt.Clean_up_partition, 'Dont do cleanup')
    matlabbatch{1}.spm.tools.oldseg.output.cleanup = 0;
elseif strcmp(opt.Clean_up_partition, 'Light Clean')
    matlabbatch{1}.spm.tools.oldseg.output.cleanup = 1;
elseif strcmp(opt.Clean_up_partition, 'Thorough Clean')
    matlabbatch{1}.spm.tools.oldseg.output.cleanup = 2;
end

%% question : tpm obligatoire? seulement Gris blanc??
if isfield(files_in, 'In2') && ~isempty(files_in.In2)
   matlabbatch{1}.spm.tools.oldseg.opts.tpm(1) =  {[files_in.In2{1}, ',1']};
end
if isfield(files_in, 'In3') && ~isempty(files_in.In3)
   matlabbatch{1}.spm.tools.oldseg.opts.tpm(2) =  {[files_in.In3{1}, ',1']};
end
if isfield(files_in, 'In4') && ~isempty(files_in.In4)
   matlabbatch{1}.spm.tools.oldseg.opts.tpm(3) =  {[files_in.In4{1}, ',1']};
end     
matlabbatch{1}.spm.tools.oldseg.opts.tpm =  matlabbatch{1}.spm.tools.oldseg.opts.tpm';
%% checker si bien rentrer (4 nombres...))
matlabbatch{1}.spm.tools.oldseg.opts.ngaus = str2num(opt.Gaussians_per_class)'; 
switch opt.Affine_regularisation
    case 'No Affine Registration'
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = '';
    case 'ICBM space template - European brains'
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'mni';
    case 'ICBM space template - East Asian brains'
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'eastern';
    case 'Average sized template'
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'subj';
    case 'No regularisation'
        matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'none';
end
matlabbatch{1}.spm.tools.oldseg.opts.warpreg = str2double(opt.Warping_regularisation);
matlabbatch{1}.spm.tools.oldseg.opts.warpco = str2double(opt.Warp_frequency_cutoff);

switch opt.Bias_regularisation
    case 'no regularisation (0)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0;
    case 'extremely light regularisation (0.00001)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.00001;
    case 'very light regularisation (0.0001)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.0001;
    case 'light regularisation (0.001)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.001;
    case 'medium regularisation (0.01)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.01;
    case 'heavy regularisation (0.1)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.1;
    case 'very heavy regularisation (1)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 1;
    case 'extremely heavy regularisation (10)'
        matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 10;
end

switch opt.Bias_FWHM
    case '30mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 30;
    case '40mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 40;
    case '50mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 50;
    case '60mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 60;
    case '70mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 70;
    case '80mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 80;
    case '90mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 90;
    case '100mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 100;
    case '110mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 110;
    case '120mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 120;
    case '130mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 130;
    case '140mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 140;
    case '150mm cutoff'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 150;
        % a vérifier !!
    case 'No correction'
        matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = Inf;
end

matlabbatch{1}.spm.tools.oldseg.opts.samp = str2double(opt.Sampling_distance);
matlabbatch{1}.spm.tools.oldseg.opts.msk = {'/Users/GIN-Eq5/Data_non_synchro/code_python/monai/data/test/C:\Faculty\TestAudrey\TemplateAndAtlas\SIGMA_ExVivo_Brain_Mask.nii,1'};
%%

if isfield(files_in, 'In5') && ~isempty(files_in.In4)
   matlabbatch{1}.spm.tools.oldseg.opts.msk =  {[files_in.In4{1}, ',1']};
else
     matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};
end 

jobs = repmat(matlabbatch, 1, 1);
inputs = cell(0, 1);
for crun = 1:1
end
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('run', jobs, inputs{:});

% rename the output file in order to match with what user expect
if isfield(files_in, 'In2') && ~isempty(files_in.In2) && ~strcmp(opt.Ouput_grey_Matter_opt, 'None')
    [path,~,ext] = fileparts(files_out.In1{1});
    movefile(fullfile(path, ['c1', char(opt.Table_in.Filename(1)),ext]), files_out.In1{1});
    copyfile(strrep(files_in.In1{1},'.nii','.json'),  strrep(files_out.In1{1},'.nii','.json') )
    %% Json processing
    if isfile(strrep(files_out.In1{1},'.nii','.json'))
        J = ReadJson(strrep(files_out.In1{1},'.nii','.json'));     
        J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
        WriteJson(J, strrep(files_out.In1{1},'.nii','.json'))
    end
end
if isfield(files_in, 'In3') && ~isempty(files_in.In3) && ~strcmp(opt.Ouput_white_Matter_opt, 'None')
   [path,~,ext] = fileparts(files_out.In2{1});
    movefile(fullfile(path, ['c2', char(opt.Table_in.Filename(1)),ext]), files_out.In2{1});
    copyfile(strrep(files_in.In1{1},'.nii','.json'),  strrep(files_out.In2{1},'.nii','.json') )
    %% Json processing
    if isfile(strrep(files_out.In2{1},'.nii','.json'))
        J = ReadJson(strrep(files_out.In2{1},'.nii','.json'));     
        J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
        WriteJson(J, strrep(files_out.In2{1},'.nii','.json'))
    end
end
if isfield(files_in, 'In4') && ~isempty(files_in.In4) && ~strcmp(opt.Ouput_csf_Matter_opt, 'None')
    [path,~,ext] = fileparts(files_out.In3{1});
    movefile(fullfile(path, ['c3', char(opt.Table_in.Filename(1)),ext]), files_out.In3{1});
    copyfile(strrep(files_in.In1{1},'.nii','.json'),  strrep(files_out.In3{1},'.nii','.json') )
    %% Json processing
    if isfile(strrep(files_out.In3{1},'.nii','.json'))
        J = ReadJson(strrep(files_out.In3{1},'.nii','.json'));     
        J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
        WriteJson(J, strrep(files_out.In3{1},'.nii','.json'))
    end
end

if ~strcmp(fullfile(path, [char(opt.Table_in.Filename(1)), '_seg_inv_sn.mat']),    files_out.In4{1})
    copyfile(fullfile(path, [char(opt.Table_in.Filename(1)), '_seg_inv_sn.mat']),  files_out.In4{1})
    copyfile(fullfile(path, [char(opt.Table_in.Filename(1)), '_seg_sn.mat']),  files_out.In5{1})
end

% add json file to .mat files
copyfile(strrep(files_in.In1{1},'.nii','.json'),  strrep(files_out.In4{1},'.mat','.json') )
%% Json processing
if isfile( strrep(files_out.In4{1},'.mat','.json'))
    J = ReadJson(strrep(files_out.In4{1},'.mat','.json'));
    J = KeepModuleHistory(J, struct('files_in', files_in, 'files_out', files_out, 'opt', opt, 'ExecutionDate', datestr(datetime('now'))), mfilename);
    WriteJson(J, strrep(files_out.In4{1},'.mat','.json'))
end
copyfile(strrep(files_out.In4{1},'.mat','.json'),  strrep(files_out.In5{1},'.mat','.json') )



%