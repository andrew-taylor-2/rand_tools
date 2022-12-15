function simple_regression(regressor,weights,outbasename,n,atlasindex,opt_string,atlasfn,full_perfusion_image_fn,cmask_fn,use_partial_data,all_masks_fn,num_subs,num_timepoints,opts)
%regressor: should be a Nx1 element numerical matrix (where N=number
%subjects). Recently added support for multiple columns, with the columns
%after the first being "nuisance variables" AKA in the model but not tested
%for effect. HOWEVER, the row/column check in the beginning of this
%function won't help if you enter these as multiple rows (isrow() only
%works on 1dim) so just make sure yourself

%weights: 1xTP element numerical matrix, weights for perfusion images
%(where TP=number time points)
%outbasename: char array, don't add a file extension.
%n: number of permutations in randomise. 5000 is standard and the default but takes a
%really long time.
%atlasfn: the atlas filename from which to grab a region
%atlasindex: the value of the area you want to extract from the atlas.
%could be a vector.
%fullperfusionimage: filename of that big perfusion image. just use
%default.

%fullperfusionimage, all_masks_fn must have all time points for subject 1
%first, then all time points for subject 2, etc.

% The functional regions could have the same name with different extents. I
% should make those names incorporate the outbasename


%try to load .mat config file
cfg_fn=fullfile(fileparts(mfilename('fullpath')),'config.mat');
if exist(cfg_fn,'file')
    config=load(cfg_fn);
    has_cfg=1;
else
%     error('you haven''t made a config.mat file with paths to data -- ')
    has_cfg=0;
end


%% sanitize inputs

if ( ~exist('full_perfusion_image_fn','var') || isempty(full_perfusion_image_fn) ) && has_cfg
    full_perfusion_image_fn=config.full_functional_image_fn;
elseif ( ~exist('full_perfusion_image_fn','var') || isempty(full_perfusion_image_fn) ) && ~has_cfg
    error('no fullperfusionimage was input and cfg file cannot be found')
%     full_perfusion_image_fn='/home/second/Desktop/new_perfusion/Data/Perfusion_all_timepoints_subs_in_MNI.nii';
end


if ( ~exist('all_masks_fn','var') || isempty(all_masks_fn) ) && has_cfg
    all_masks_fn=config.full_data_coverage_image_fn;
elseif ( ~exist('all_masks_fn','var') || isempty(all_masks_fn) ) && ~has_cfg && use_partial_data
    error('no all masks fn was input and cfg file cannot be found')
%     all_masks_fn='/home/second/Desktop/new_perfusion/Data/66masks.nii';
end

if ( ~exist('atlasfn','var') || isempty(atlasfn) ) && has_cfg && ( isfield(config,'fsldir') && ~isempty(config.fsldir) )
    fsldir=config.fsldir;
    atlasfn=fullfile(fsldir,'data','atlases','Juelich','Juelich-maxprob-thr25-2mm.nii.gz');
    
elseif ( ~exist('atlasfn','var') || isempty(atlasfn) )
    fsldir=getenv('FSLDIR');
    if ~isempty(fsldir)
        atlasfn=fullfile(fsldir,'data','atlases','Juelich','Juelich-maxprob-thr25-2mm.nii.gz');
    else
%         error(' couldn''t find your atlas_fn')
        atlasfn='/usr/local/fsl/data/atlases/Juelich/Juelich-maxprob-thr25-2mm.nii.gz';
    end
end
assert(logical(exist(atlasfn,'file')))

if ~exist('opt_string','var') || isempty(opt_string)
    opt_string='';
end
if ~exist('n','var') || isempty(n)
    n='5000';
elseif isscalar(n)
    n=num2str(n);
end

if ( ~exist('cmask_fn','var') || isempty(cmask_fn) ) && has_cfg
    cmask_fn=config.consensus_mask_fn;
elseif ( ~exist('cmask_fn','var') || isempty(cmask_fn) ) && ~has_cfg
    %not gonna give an error here, because this variable might not be
    %needed if use_partial_data

%     cmask_fn='/home/second/Desktop/new_perfusion/Data/common_mask_all_in_MNI.nii.gz';
end

if ~isrow(weights)
    weights=weights';
end

if isrow(regressor)
    regressor=regressor';
end

if ~exist('use_partial_data','var') || isempty(use_partial_data)
    use_partial_data=0;
end

if ~exist('num_subs','var') || isempty(num_subs)
    num_subs=11;
end


if ~exist('num_timepoints','var') || isempty(num_timepoints)
    num_timepoints=6;
end

%tfce option
if ~exist('opts','var') || ~isfield(opts,'t')
    opts.t=1;
end

[pth,nme,~]=fileparts(outbasename);

%% make the design and contrast matrices
demean=@(x) x-mean(x(:));

%separate nuisance and regressor
if size(regressor,2)>1
    disp('simple_regression thinks that you are including nuisance variables in your design -- if you''re not, consult with andrew lol')
    nuisance=regressor(:,2:end);
    regressor=regressor(:,1);
    
    for i=1:size(nuisance,2)
        nuisance(:,i)=demean(nuisance(:,i));
    end
    
    numnuis=size(nuisance,2);
    zs=zeros(2,numnuis);
        
    des=[ones(length(regressor),1),demean(regressor),nuisance];
    con=[[0 1; 0 -1] zs];
    
else %if you just have the one covariate
    des=[ones(length(regressor),1),demean(regressor)];
    con=[0 1; 0 -1];
    
end


%if output folder isnt' a folder yet, make it
mkdir(pth)

%write text file
[matfn,confn]=make_design(des,con,[outbasename '_mat.txt'],[outbasename '_con.txt']);

randimgfn=perf_model(full_perfusion_image_fn,weights,pth,[nme '_perfusion_image'],num_subs,num_timepoints);

% isolate binary region from atlas
atlas=d2n2s(atlasfn);
functionalarea=just_these(atlas,atlasindex,1); % recently added a binarize option to just_these --
% I'm not sure if randomise treats masks differently based on their values
% or if it just binarizes based on ~=0 itself. So it's safest to binarize.

%come up with a good name for the region
[~,atlname,~]=fileparts(atlasfn);
if numel(atlname)>6
    atlseg=atlname(1:7);
else
    atlseg=atlname;
end
atlseg=[atlseg strrep(num2str(atlasindex),' ','_')];

%I think find_the_best_mask() should come here.
if use_partial_data
    
    %get the various masks
    [bin_mask,lmask_for_randomise,num_subs_image]=find_the_best_mask(functionalarea,weights,all_masks_fn,num_subs,num_timepoints,num_subs-4);
    
    %then turn your .mats into new .mats with setup_masks
    [outmat,outcon,rand_cmd_addition,outmsg]=setup_masks(matfn,confn,outbasename,lmask_for_randomise);
    
    %to avoid confusion
    delete(matfn)
    delete(confn)
    
    matfn=outmat;
    confn=outcon;
    
    %write bin_mask and num_subs_image for 1 -m option randomise and 2. qc.
    ROI_fn=d2n2s_write(bin_mask,pth,[nme '_' atlseg '_partial_masked'],'dt',[0 2],'del',1)
    d2n2s_write(num_subs_image,pth,[nme '_num_subs_in_test'],'dt',[0 2],'del',1)

    
else %actually, i think the rest of the script should end up in this conditional
    rand_cmd_addition='';
    
    
    %load the FOV inclusivity image, which we'll refer to as a "common mask"
    cmask=d2n2s(cmask_fn);
    
    %now we want to set our motor area to 0 anywhere that our common mask has
    %1 (after the inversion we just did). On the next line is an assignment that only assigns to array elements
    %for which logical_mask has a true value.
    functionalarea.img(~logical(cmask.img))=0;
    
    if isequal(unique(functionalarea.img(:)),0)
        error('the image for the region being tested in randomise has no non-zero values -- likely the entire region was outside the common mask')
    end
    
    
    %write the new image and assign the filename to a variable
    ROI_fn=d2n2s_write(functionalarea,pth,atlseg,'dt',[0 2],'del',1);
    %note that there really isn't much left of this area after we mask it --
    %maybe we can't do any tests in this area.
    
end


%make a new output base name

%tfce or no
if ~opts.t
    tfce_string='';
else
    tfce_string=' -T';
end

system(['randomise -i ' randimgfn ' -o ' outbasename ' -d ' matfn ' -t ' confn ' -m ' ROI_fn ' --uncorrp ' opt_string rand_cmd_addition ' -n ' n tfce_string ' > ' outbasename '_log.txt 2>' outbasename '_err.txt &'])
disp(sprintf(['\n **randomise is running in the background** \n **and its output is going to ' outbasename '_log.txt**']))
    


% %also gzip images lol
%
%
% %load the FOV inclusivity image, which we'll refer to as a "common mask"
% cmask=d2n2s(cmask_fn);
%
% %now we want to set our motor area to 0 anywhere that our common mask has
% %1 (after the inversion we just did). On the next line is an assignment that only assigns to array elements
% %for which logical_mask has a true value.
% motorarea.img(~logical(cmask.img))=0;
%
% if isequal(unique(motorarea.img(:)),0)
%     error('the image for the region being tested in randomise has no non-zero values -- likely the entire region was outside the common mask')
% end
%
% %come up with a good name for the region
% [~,atlname,~]=fileparts(atlasfn);
% if numel(atlname)>6
%     atlseg=atlname(1:7);
% else
%     atlseg=atlname;
% end
% atlseg=[atlseg strrep(num2str(atlasindex),' ','_')];
%
% %write the new image and assign the filename to a variable
% marea_fn=d2n2s_write(motorarea,pth,atlseg,'dt',[0 2],'del',1);
% %note that there really isn't much left of this area after we mask it --
% %maybe we can't do any tests in this area.
%
% %make a new output base name
% system(['randomise -i ' randimgfn ' -o ' outbasename ' -d ' matfn ' -t ' confn ' -m ' marea_fn ' --uncorrp ' opt_string rand_cmd_addition ' -n ' n ' -T > ' outbasename '_log.txt &'])
% disp(sprintf(['\n **randomise is running in the background** \n **and its output is going to ' outbasename '_log.txt**']))
