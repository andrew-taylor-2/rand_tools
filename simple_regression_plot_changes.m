function varargout=simple_regression_plot_changes(regressor,weights,outbasename,~,~,~,~,full_perfusion_image_fn,~,use_partial_data,all_masks_fn,num_subs,num_timepoints,sig_region_fn,perf_label,beh_label,sans_vec,graphcase,opts)
%regressor: should be a Nx1 element numerical matrix (where N=number
%subjects)
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

do=get_anonymous_functions;

if ~exist('graphcase','var')
    graphcase='';
end

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

if ( ~exist('fullperfusionimage','var') || isempty(full_perfusion_image_fn) ) && has_cfg
    full_perfusion_image_fn=config.full_functional_image_fn;
elseif ( ~exist('fullperfusionimage','var') || isempty(full_perfusion_image_fn) ) && ~has_cfg
    error('no fullperfusionimage was input and cfg file cannot be found')
%     full_perfusion_image_fn='/home/second/Desktop/new_perfusion/Data/Perfusion_all_timepoints_subs_in_MNI.nii';
end


if ( ~exist('all_masks_fn','var') || isempty(all_masks_fn) ) && has_cfg
    all_masks_fn=config.full_data_coverage_image_fn;
elseif ( ~exist('all_masks_fn','var') || isempty(all_masks_fn) ) && ~has_cfg
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
% assert(exist(atlasfn,'file'))
% 
% if ~exist('opt_string','var') || isempty(opt_string)
%     opt_string='';
% end
% if ~exist('n','var') || isempty(n)
%     n='5000';
% elseif isscalar(n)
%     n=num2str(n);
% end
% 
% if ( ~exist('cmask_fn','var') || isempty(cmask_fn) ) && has_cfg
%     cmask_fn=config.consensus_mask_fn;
% elseif ( ~exist('cmask_fn','var') || isempty(cmask_fn) ) && ~has_cfg
%     %not gonna give an error here, because this variable might not be
%     %needed if use_partial_data
% 
% %     cmask_fn='/home/second/Desktop/new_perfusion/Data/common_mask_all_in_MNI.nii.gz';
% end

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

if ~exist('opts','var') || ~isfield(opts,'trend')
    opts.trend=0;
end

if ~exist('opts','var') || ~isfield(opts,'thresh')
    opts.thresh=.949999;
end

[pth,nme,~]=fileparts(outbasename);

%% make the design and contrast matrices
% demean=@(x) x-mean(x(:));
% 
% des=[ones(length(regressor),1),demean(regressor)];
% 
% con=[0 1; 0 -1];

%if output folder isnt' a folder yet, make it
% mkdir(pth)

%write text file
% [matfn,confn,~,~]=make_design(des,con,[outbasename '_mat.txt'],[outbasename '_con.txt']);

regressor=regressor(:,1);
[perf_rois,roi,inds]=perf_model_region_time_course(sig_region_fn,full_perfusion_image_fn,weights,num_subs,num_timepoints,opts.thresh);
%indexing doesn't have to be crazy bc there's a mask for each sub, not for
%each time point. Or maybe there is a mask for each TP....... no no just
%kidding there isn't 

%also I won't need this atlas code bc the region is necessarily within the
%atlas
% 
% % isolate binary region from atlas
% atlas=d2n2s(atlasfn);
% functionalarea=just_these(atlas,atlasindex,1); % recently added a binarize option to just_these --
% % I'm not sure if randomise treats masks differently based on their values
% % or if it just binarizes based on ~=0 itself. So it's safest to binarize.
% 
% %come up with a good name for the region
% [~,atlname,~]=fileparts(atlasfn);
% if numel(atlname)>6
%     atlseg=atlname(1:7);
% else
%     atlseg=atlname;
% end
% atlseg=[atlseg strrep(num2str(atlasindex),' ','_')];

%I think find_the_best_mask() should come here.
if use_partial_data
    
    %make a functional area that includes everything.
    functionalarea=perf_rois{1}{1}; %this perf_rois should always exist....
    functionalarea.img=double(roi); 
    
    %get the various masks
    [bin_mask,lmask_for_randomise,num_subs_image]=find_the_best_mask(functionalarea,weights,all_masks_fn,num_subs,num_timepoints,num_subs-4);
    %functionalarea used to include areas where the mask was, but we
    %shouldn't need that bc the whole region is necessarily in it. 
    
    
    %then turn your .mats into new .mats with setup_masks
%     [outmat,outcon,rand_cmd_addition,outmsg]=setup_masks(matfn,confn,outbasename,lmask_for_randomise);
%we don't need that line at all in this script    

    %to avoid confusion
%     delete(matfn)
%     delete(confn)
    
%     matfn=outmat;
%     confn=outcon;
    
    %write bin_mask and num_subs_image for 1 -m option randomise and 2. qc.
%     ROI_fn=d2n2s_write(bin_mask,pth,[nme '_' atlseg '_partial_masked'],'dt',[0 2],'del',1)
%     d2n2s_write(num_subs_image,pth,[nme '_num_subs_in_test'],'dt',[0 2],'del',1)

%don't need to write for this script.

    
else %actually, i think the rest of the script should end up in this conditional
    rand_cmd_addition='';
    
    functionalarea=d2n2s(sig_region_fn,'no','bvecbvaljson');
    functionalarea.img(functionalarea.img<opts.thresh)=0;
    
    
    %write the new image and assign the filename to a variable
%     ROI_fn=d2n2s_write(functionalarea,pth,atlseg,'dt',[0 2],'del',1);
%don't need to write for this script

end

%okay so normally randomise does masking for us, so I'm gonna have to write
%a new function for masking now. 
%   we'll use either functionalarea (no partial) or bin_mask and lmask.
%   gonna have to use inds.(blah) to figure out what of lmask to use right?
% is lmask one per subject or one per relevant timepoint?

%mask EACH perf_roi and then get average for EACH element of perf_roi

%i have these index changy things
% inds.stp_2_1dim{i}{tp}=(i-1)*num_timepoints+(j);
% inds.stp_2_ij{i}{tp}=[i,j];


%sub-mask the perf roi:
num_tps=sum(logical(weights));
for i=1:num_subs
    
    %invert lmask
    lmask_for_randomise(i).img(:)=~lmask_for_randomise(i).img; %stephen cobeldick says it's cool
        
    %find an AND of lmask and sig region/functionalarea
    lmask_for_randomise(i).img=lmask_for_randomise(i).img & functionalarea.img; % no need to repmat since we're doing it for every sub    repmat(functionalarea.img,[1 1 1 num_subs])
    
    for tp=1:num_tps
        %for each tp find the mean
        data{i}{tp}=perf_rois{i}{tp}.img(logical(lmask_for_randomise(i).img));
        means{i}(tp)=mean(data{i}{tp},'all');
    end
end

%% begin the graphing





xaxlabels={'BDC13','BDC7','HDT7','HDT29','R5','R12'};
xaxlabels(~logical(weights))=[];
%okay now we already have what to label the stuff with. 
cvecs=distinguishable_colors(num_subs);


%% graphing cases

switch graphcase
    case {'','perf-beh-correlation'}
        
        %% add brain 3 views
        
        figure;
        mnit1=d2n2s(fullfile(fsldir,'data','standard','MNI152_T1_2mm.nii.gz'));
        regg=d2n2s(sig_region_fn,'no','bvecbvaljson');
        
        [~,maxind]=max(regg.img(:));
        [x,y,z]=ind2sub(size(regg.img),maxind);
        
        %thresh regg
        regg.img(regg.img<opts.thresh)=0;
        
        %subplot
        t=tiledlayout(2,4,'TileSpacing','compact')
        s=@squeeze;
        %z
        nexttile;
        imout=imoverlay(mat2gray(mnit1.img(:,:,z)),regg.img(:,:,z),'yellow');
        imshow(imout);
        axis square
        
        %y
        nexttile;
        imout=imoverlay(s(mat2gray(mnit1.img(:,y,:))),s(regg.img(:,y,:)),'yellow');
        imshow(imout);
        axis square
        
        %x
        nexttile(5);
        imout=imoverlay(s(mat2gray(mnit1.img(x,:,:))),s(regg.img(x,:,:)),'yellow');
        imshow(imout);
        axis square
        
        
        
        %% add contrast score
        %need to apply the contrasts to these
        
        nexttile(3,[2 2])
        
        
        relweights=weights;
        relweights(weights==0)=[];
        for i=1:num_subs
            summ{i}=zeros(size(data{i}{tp}));
            for tp=1:num_tps
                summ{i}=summ{i}+data{i}{tp}.*relweights(tp);
            end
            conmeans(i)=mean(summ{i},'all');
        end
        
        
        %% plot
        hold on
        for i=1:length(regressor)
            if sans_vec(i)==1
                scatter(conmeans(i),regressor(i),50,cvecs(i,:),'Marker','+')
            else
                scatter(conmeans(i),regressor(i),50,cvecs(i,:),'Marker','o')
            end
        end
        hold off
        
        ax3=gca
        set(get(ax3,'XLabel'),'String',['Perfusion Contrast Value: ' perf_label])
        set(get(ax3,'YLabel'),'String',['Behavior Contrast Value: ' beh_label])
        
        % add trends?
        if opts.trend==2
            mdlSANS=fitlm(conmeans(sans_vec==1),regressor(sans_vec==1));
            mdlNSANS=fitlm(conmeans(sans_vec~=1),regressor(sans_vec~=1));
            
            
            refline(mdlSANS.Coefficients.Estimate(2),mdlSANS.Coefficients.Estimate(1))
            refline(mdlNSANS.Coefficients.Estimate(2),mdlNSANS.Coefficients.Estimate(1))
            
        elseif opts.trend==1
            mdl=fitlm(conmeans,regressor);
            coefs=mdl.Coefficients.Estimate;
            refline(coefs(2),coefs(1))
        else
        end
        
        
        
        
        
    case 'perf-timecourse'
        
        %just graph "data" by time point
        hold on
        for i=1:num_subs
            tpdatameans{i}=cellfun(@mean,data{i});
            if sans_vec(i)==1
                plot(tpdatameans{i},'Color',cvecs(i,:),'LineStyle','--')  %using plot(y) syntax
            else
                plot(tpdatameans{i},'Color',cvecs(i,:),'LineStyle','-')  %using plot(y) syntax
            end
        end
        xticklabels(xaxlabels)
        
        %out
        varargout{1}=tpdatameans;
               
        
    case 'dontgraph_returndata'
        
        for i=1:num_subs
            tpdatameans{i}=cellfun(@mean,data{i});
        end
        
        varargout{1}=tpdatameans;
        varargout{2}=cvecs;
        varargout{3}=sans_vec;
                
        
end



% set(ax3,'XTick',[])

%ADD LABEL FOR SANS SUBS:
% I should use different markers 
%eg 'Marker','+' OR 'Marker','o' for sans,
%nonsans

%could i (instead of using whole significant region)
% just use the most significant voxel? (this is too sophisticated to
% implement now but also what if i used the biggest region where the most
% subs had data)

%% okay now add the part where we display the region on MNI
% 
% % COPY AND PASTED SPM GEMS FOLLOWS
% % Make sure to first clear the graphics window
%  
% % Select images
% % Pbck = spm_get(1,'*.img','Select background image')
% Pbck=do.gunzip_and_rename(do.copy_and_rename(fullfile(fsldir,'data','standard','MNI152_T1_2mm.nii.gz'),'mnit12mm.nii.gz'));
% % Psta = spm_get(1,'*.img','Select statistic image')
% Psta=do.gunzip_and_rename_no_delete(sig_region_fn);
%  
% % Set the threshold value
% Th   = 0.95;  
%  
% % Create a new image where all voxels below Th have value NaN
% PstaTh = do.append(Psta,'nan4disp');
% spm_imcalc(Psta,PstaTh,'i1+(0./(i1>=Th))',{[],[],spm_type('float')},Th);
%  
% % Display!
% clear global st
% ff=figure;
% spm_orthviews('image',Pbck,[0.05 0.05 0.9 0.9],ff); %reposition using the line that happens later
% spm_orthviews('addimage',1,PstaTh)
%  
% % Possibly, set the crosshairs to your favorite location
% 
% %find most significant voxel, reposition to it
% %cant' believe im loading this again
% regg=d2n2s(sig_region_fn);
% [~,maxind]=max(regg.img(:));
% [x,y,z]=ind2sub(size(regg.img),maxind)
% ijk=regg.hdr.mat*[x,y,z,1]';
% clear global st
% spm_orthviews('reposition',[ijk(1) ijk(2) ijk(3)])
% delete(Pbck)
% % END GEMS 








% mnit1=d2n2s(fullfile(fsldir,'data','standard','MNI152_T1_2mm.nii.gz'));

% 
% function make_panel(t1,axes1,backgd,foregd,pos)
% 
% im = imagesc(axes1,backgd); 
% im.AlphaData = 0.5; % change this value to change the background image transparency 
% axis square; 
% hold all; 
% %plot second data 
% axe2 = nexttile(pos); %okay, this line isn't even showing the second image
% % set(axe2,'Parent',t1) % this command sticks the axes in the first position of the tiled layout
% % %how could i set it to match the position of the other axes;
% 
% im1 = imagesc(axe2,foregd); 
% im1.AlphaData = 0.5; % change this value to change the foreground image transparency 
% axis square; 
% %link axes 
% linkaxes([axes1,axe2]) 
% %%Hide the top axes 
% axe2.Visible = 'off'; 
% axe2.XTick = []; 
% axe2.YTick = []; 
% %add differenct colormap to different data if you wish 
% colormap(axes1,'gray') 
% colormap(axe2,'jet') 
% %set the axes and colorbar position 


% 
% if false
%     % tiled layout ideas
%     tiledlayout(1,5)
%     nexttile(1,[1 3])
%     
%     %% this block is specific to the line graph idea
%     %get the distinguishable colors here
%     for i=1:num_subs
%         
%         plot(1:num_tps,means{i},'Color',cvecs(i,:),'LineWidth',3.5)
%         hold on
%     end
%     hold off
%     
%     ax=gca
%     set(ax,'XTick',1:num_tps)
%     set(ax,'XTickLabel',xaxlabels)
%     set(ax,'XLim',[1-.3 num_tps+.3])
%     
%     
%     %% okay, now add the stuff for beh scores
%     nexttile
%     
%     scatter(ones(1,length(regressor)),regressor,35,cvecs,'filled')
%     ax2=gca
%     set(ax2,'XTick',[])
% end