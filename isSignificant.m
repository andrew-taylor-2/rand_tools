function proportion=isSignificant(image_fn)
% if isstruct(image_fn)
%     image=image_fn;
% else
    image=d2n2s(image_fn,'no','bvecbvaljson');
% end
proportion=sum(image.img(:)>=.95)/numel(image.img);
end