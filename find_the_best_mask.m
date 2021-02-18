function find_the_best_mask(weights_raw,mask_img_all)

%tasks: make a best mask 


weights_raw=[.5 .5 0 -1 0 0];


weights=repmat(weights_raw,[1 11])
% This version of weights is for a subject -major image. That is, the first
% couple images are all from the same subject. Like img_all.

% W A I T,
% I think I've been getting the meaning of "-major" wrong. BUt whatever.
% context clues.

for i=1:11
    for j=1:6
        disp((i-1)*6+(j))
        i
        j
    end
end

%haven't run this next part yet
clear wimg_all sub_sum %clearing sub sum is dumb given the next line but whatever
sub_sum=zeros(size(img_all(1).img,1),size(img_all(1).img,2),size(img_all(1).img,3),11);
for i=1:11 %subs
    %sub_sum(:,:,:,i)=zeros(size(img_all(1).img,1),size(img_all(1).img,2),size(img_all(1).img,3));
    for j=1:6 %time points 
        %random labeling line to help keep things together.
        %img_all((i-1)*6+(j)).label=index_path_from_end(perfglob((i-1)*6+(j)).folder,2:3);
        
        wimg_all((i-1)*6+(j))=img_all((i-1)*6+(j)); %make them equal for convenience, then change the image
        wimg_all((i-1)*6+(j)).img=img_all((i-1)*6+(j)).img.*weights_raw(j);
        
        sub_sum(:,:,:,i)=sub_sum(:,:,:,i)+wimg_all((i-1)*6+(j)).img;
        %sub_sum_label{i,j}=wimg_all((i-1)*6+(j)).label;
    end
end
% 
% for i=1:11
%     for j=1:6
%         disp((i-1)*6+(j))
%         disp(i)
%         disp(j)
%         sub_sum_label{i,j}
%     end
% end

for i=1:11
temp(i)=img_all(1);
temp(i).img=sub_sum(:,:,:,i);
end
inimg27=d2n2s_write(temp,'/Users/dr.robertslab/Andrew/Projects/NASA/d_redo/final_base_to_hdt29','base_to_hdt29_drop')
