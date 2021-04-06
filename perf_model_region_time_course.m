function perf_rois=perf_model_region_time_course(region_fn,in_fn,weights_raw,num_subs,num_timepoints)

%I can't just get the mean because I don't know how it will be
%partial masked.......

img_all=d2n2s(in_fn);
region=d2n2s(region_fn);
roi=logical(region.img);

%make templates
region.img=zeros(size(region.img));
temp=repmat(region,[num_subs.*num_timepoints 1]);

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
            theregion=img_all((i-1)*num_timepoints+(j)).img(roi);
            
            perf_rois{i}(tp)=mean(theregion,'all');
        end
        
    end
end

