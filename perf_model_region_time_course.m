function [perf_rois,roi,inds]=perf_model_region_time_course(region_fn,in_fn,weights_raw,num_subs,num_timepoints)
do=get_anonymous_functions;
%what do i want out of this? I think it's gonna be just the relevant area
%from the region out of img_all. but it's got to be a full image/nii bc i
%with later mask it some more

%I can't just get the mean because I don't know how it will be
%partial masked.......

img_all=d2n2s(in_fn);
region=d2n2s(region_fn);
roi=region.img>0.94999;

%make templates
region.img=zeros(size(region.img));
[perf_rois{1}{1:num_timepoints}]=deal(do.curly(repmat({region},[num_timepoints 1]),':'));
perf_rois=repmat(perf_rois,[num_subs 1]);
%perf rois is now the right size and everything
%perf_rois{num subs}{time points}

%get the right indices
weights_bool=logical(weights_raw);

weights=repmat(weights_bool,[1 num_subs]); %am i going to use this anymore?

%loop uses i and j so let's make an array with i and j
%i can use i directly but instead of j i want what position the weight is
%in. I can just use a counter i guess.


%% new loop
for i=1:num_subs %subs
    tp=0; %reset so that the first match is 1 again 
    
    for j=1:num_timepoints %time points
        
        
        if weights_bool(j) %it's a tp we're looking at
            tp=tp+1;
%             theregion=img_all((i-1)*num_timepoints+(j)).img(roi);
            img_all((i-1)*num_timepoints+(j)).img(~roi)=0;
%             perf_rois{i}{tp}=mean(theregion,'all');
            perf_rois{i}{tp}.img=img_all((i-1)*num_timepoints+(j)).img;
            
            inds.stp_2_1dim{i}{tp}=(i-1)*num_timepoints+(j);
            inds.stp_2_ij{i}{tp}=[i,j];
            
        end
        
    end
end

%I'm guessing I'm going to have to do that tp indexing later, too? 
