function [bin_mask,lmask_for_randomise,num_subs_image]=find_the_best_mask(motorarea,weights,all_masks_fn)
% gets you randomise mask inputs and images for QC

%load all_masks_fn
all_masks=d2n2s(all_masks_fn,'no','bvalbvecjson');

%grab one obj as a template
template=all_masks(1);

%sanitize input
weights=logical(weights);

%% sum across timepoints

sub_sum=zeros(size(all_masks(1).img,1),size(all_masks(1).img,2),size(all_masks(1).img,3),11);
for i=1:11 %subs
    for j=1:6 %time points 
        
        %make them equal for convenience
        wall_masks((i-1)*6+(j))=all_masks((i-1)*6+(j)); 
        
        %here, stuff will only be counted if weights is not 0/false
        wall_masks((i-1)*6+(j)).img=all_masks((i-1)*6+(j)).img.*weights_raw(j);
        
        %get a sum across timepoints for each subject
        sub_sum(:,:,:,i)=sub_sum(:,:,:,i)+wall_masks((i-1)*6+(j)).img;
        
    end
end


%% make the output images

%FIRST: we must make a within-subject consensus mask. we can't have partial
%time-point data, because it would be weighted super wrong.

%how many timepoints is all the relevant timepoints?
consensus_thresh=sum(logical(weights(1:6)));

%make a XxYxZx11 boolean image (the vox value isn't important bc it's
%just the number of relevant timepoints or 0)
sub_sum=binarize_to_value(sub_sum,1,consensus_thresh);

%get sum of subjects with full data in voxel
num_subs_image=sum(sub_sum,4);

%lets get this in just the ROI, then 0 any values under 7. That's one
%output. Then we use that as a boolean to cut out values from the 4D image.

%get in just the ROI. 
num_subs_image(~logical(motorarea.img))=0;

%0 any values under 7
num_subs_image(num_subs_image<7)=0;
%that's an output

%binarize to see full region tested
bin_mask=double(logical(num_subs_image));
%that's an output

%get the 4d image which shows where data is available (0 value) for each
%subject
lmask_for_randomise=double(~logical(sub_sum(logical(num_subs_image))));
%that's an output

%% put images in objects

%metadata is same for all so assign template
[num_subs_obj,bin_mask_obj,lmask_obj]=deal(template);

%assign correct voxel data
num_subs_obj.img=num_subs_image;
bin_mask_obj.img=bin_mask;
lmask.img=lmask_for_randomise;