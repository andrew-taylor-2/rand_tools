function randimgfn=perf_model_region_time_course(region_fn,in_fn,weights_raw,folderr,filenamee,num_subs,num_timepoints)
% img_all=d2n2s(dir('/home/second/Desktop/new_perfusion/Data/Perf*.nii'));

img_all=d2n2s(in_fn);
region=d2n2s(region_fn);

%get the right indices
weights_bool=logical(weights_raw)

weights=repmat(weights_bool,[1 num_subs]);

%loop uses i and j so let's make an array with i and j
%i can use i directly but instead of j i want what position the weight is
%in. I can just use a counter i guess.
tp=0;

if weights_bool(j) %it's a tp we're looking at 
    tp=tp+1;
    rel_roi{i}{tp}=mean();
end


% clear wimg_all sub_sum 
% sub_sum=zeros(size(img_all(1).img,1),size(img_all(1).img,2),size(img_all(1).img,3),num_subs);
for i=1:num_subs %subs
    for j=1:num_timepoints %time points 

        wimg_all((i-1)*num_timepoints+(j))=img_all((i-1)*num_timepoints+(j)); %make them equal for convenience, then change the image
        wimg_all((i-1)*num_timepoints+(j)).img=img_all((i-1)*num_timepoints+(j)).img.*weights_raw(j);
        
        sub_sum(:,:,:,i)=sub_sum(:,:,:,i)+wimg_all((i-1)*num_timepoints+(j)).img;
    end
end

%what should the output look like?

temp=img_all(1);
temp=repmat(temp,[])

%output doesn't need to be object, I can just get means staight from here


for i=1:num_subs
temp(i)=img_all(1);
temp(i).img=sub_sum(:,:,:,i);
end
randimgfn=d2n2s_write(temp,folderr,filenamee,'del',1)
