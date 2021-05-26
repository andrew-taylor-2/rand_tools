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

if ~exist('opts','var') || ~isfield(opts,'flip')
    opts.flip=true;
end

%kludge
if ~exist('opts','var') || ~isfield(opts,'subs')
    opts.sv={'A','B','C','D','E','G','H','J','K','L','M'};
end

%take out subs if the user says to
if ~exist('opts','var') || ~isfield(opts,'exclude')
    opts.exclude=[];
end

%remove
opts.sv(opts.exclude)=[];

%see if it matches with the number of subs we've given
if numel(opts.sv)~=num_subs
    error('the number of subjects doesn''t match the vector of subject labels. \nUse opts.exclude=[excluded subjects] to indicated which subs to take out')
end



[pth,nme,~]=fileparts(outbasename);

%% make the design and contrast matrices

regressor=regressor(:,1);
[perf_rois,roi,inds]=perf_model_region_time_course(sig_region_fn,full_perfusion_image_fn,weights,num_subs,num_timepoints,opts.thresh);
%indexing doesn't have to be crazy bc there's a mask for each sub, not for
%each time point. Or maybe there is a mask for each TP....... no no just
%kidding there isn't 


%I think find_the_best_mask() should come here.
if use_partial_data
    
    %make a functional area that includes everything.
    functionalarea=perf_rois{1}{1}; %this perf_rois should always exist....
    functionalarea.img=double(roi); 
    
    %get the various masks
    [bin_mask,lmask_for_randomise,num_subs_image]=find_the_best_mask(functionalarea,weights,all_masks_fn,num_subs,num_timepoints,num_subs-4);
    %functionalarea used to include areas where the mask was, but we
    %shouldn't need that bc the whole region is necessarily in it. 
    
    
    
else %actually, i think the rest of the script should end up in this conditional
    rand_cmd_addition='';
    
    functionalarea=d2n2s(sig_region_fn,'no','bvecbvaljson');
    functionalarea.img(functionalarea.img<opts.thresh)=0;
    

end



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
        t=tiledlayout(2,13,'TileSpacing','compact')
        s=@squeeze;
        %z
        nexttile(1,[1 3]);
        imout=imoverlay(flip(mat2gray(mnit1.img(:,:,z))'),flip(regg.img(:,:,z)'),'yellow');
        imshow(imout);
        axis square
        %title('Neurological');
        
        
        %y
        nexttile(4,[1 3]);
        imout=imoverlay(flip(s(mat2gray(mnit1.img(:,y,:)))'),flip(s(regg.img(:,y,:))'),'yellow');
        imshow(imout);
        axis square
        
        %x
        nexttile(14,[1 3]);
        imout=imoverlay(flip(s(mat2gray(mnit1.img(x,:,:)))'),flip(s(regg.img(x,:,:))'),'yellow');
        imshow(imout);
        axis square
        
        %text
%         nexttile(6);
%         text(1,1,{'result: positive correlation','x data: behavior stuff'})
%         axis square
        
%         if ~isempty(opts.table)
%             opts.table
%             Result='positive correlation'
%             vox_data='perfusion drop'
%             
%             uitable('Data',Table(Result,vox_data),'ColumnName',T.Properties.VariableNames,...
%     'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%         
%         
        
        %% add contrast score
        %need to apply the contrasts to these
        
        nexttile(8,[2 6])
        
        
        relweights=weights;
        relweights(weights==0)=[];
        for i=1:num_subs
            summ{i}=zeros(size(data{i}{tp}));
            for tp=1:num_tps
                summ{i}=summ{i}+data{i}{tp}.*relweights(tp);
            end
            conmeans(i)=mean(summ{i},'all');
            
                
        end
        if opts.flip
            conmeans=-1*conmeans;
        end
        
        %% plot
        hold on
        for i=1:length(regressor)
            if sans_vec(i)==1
                h(i)=scatter(conmeans(i),regressor(i),75,cvecs(i,:),'^','filled');
            else
                h(i)=scatter(conmeans(i),regressor(i),75,cvecs(i,:),'o','filled');
            end
        end
        
        % now we make black markers for sans, non sans
        L(1)=scatter(0,0,75,[0 0 0 ],'^','filled');
        L(2)=scatter(0,0,75,[0 0 0 ],'o','filled');

        
        hold off
        
        ax3=gca
        set(get(ax3,'XLabel'),'String',[ perf_label],'FontSize',14)
        set(get(ax3,'YLabel'),'String',[ beh_label],'FontSize',14)
        
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
        
  %JOE: Idea to add a lengend but still needs some editing
  
        %legend    
        legend([h L]',[opts.sv {'SANS','nonSANS'}]','AutoUpdate','Off');
            %Gets the square/circle right but wonder if there is a way to
            %make them plain
            
         %Adding dotted lines at x=0 and y=0
         xline(0,':');
         yline(0,':');

        %now make the 0 points go away
        set(L,'Visible','off')
        
        
        
        
    case 'perf-timecourse'
        
        if opts.flip
            for i=1:num_subs
                tpdatameans{i}=-1*tpdatameans{i};
            end
        end
        
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
        
        if opts.flip
            for i=1:num_subs
                tpdatameans{i}=-1*tpdatameans{i};
            end
        end
        
        for i=1:num_subs
            tpdatameans{i}=cellfun(@mean,data{i});
        end
        
        varargout{1}=tpdatameans;
        varargout{2}=cvecs;
        varargout{3}=sans_vec;
                
        
end

