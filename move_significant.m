function move_significant(inbasefolder,test_str,outfolder,outbasefolder)

if ~exist('outbasefolder','var')
    outbasefolder='voxelwise_36';
end

% test_str='*JHUdiv*corrp*'
destination=[ outbasefolder filesep outfolder filesep 'tests_with_significant_voxels/'];
mkdir(destination)

count=0;
tstatdir=clean_dir(dir(strrep([ inbasefolder filesep test_str '*corrp*.nii.gz'],'**corrp','*corrp')));
for i=1:length(tstatdir)
    imgg=d2n2s(tstatdir(i));
    if any(imgg.img(:)>=.95)
        copyfile(fnify2(tstatdir(i)),destination)
        count=count+1;
    else
    end
end
disp('number of significant things')
disp(count)
