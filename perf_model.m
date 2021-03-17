function randimgfn=perf_model(in_fn,in_weights,folderr,filenamee,num_subs,num_timepoints)
% img_all=d2n2s(dir('/home/second/Desktop/new_perfusion/Data/Perf*.nii'));
img_all=d2n2s(in_fn);
% weights_raw=[.5 .5 -1 0 0 0];
weights_raw=in_weights;


weights=repmat(weights_raw,[1 num_subs]);
%the first couple images are all from the same subject. Like img_all.


clear wimg_all sub_sum 
sub_sum=zeros(size(img_all(1).img,1),size(img_all(1).img,2),size(img_all(1).img,3),num_subs);
for i=1:num_subs %subs
    for j=1:num_timepoints %time points 

        wimg_all((i-1)*num_timepoints+(j))=img_all((i-1)*num_timepoints+(j)); %make them equal for convenience, then change the image
        wimg_all((i-1)*num_timepoints+(j)).img=img_all((i-1)*num_timepoints+(j)).img.*weights_raw(j);
        
        sub_sum(:,:,:,i)=sub_sum(:,:,:,i)+wimg_all((i-1)*num_timepoints+(j)).img;
    end
end

for i=1:num_subs
temp(i)=img_all(1);
temp(i).img=sub_sum(:,:,:,i);
end
randimgfn=d2n2s_write(temp,folderr,filenamee,'del',1)
