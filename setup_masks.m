function [outmat,outcon,randomise_cmd_addendum,outmsg]=setup_masks(inmat,incon,outbasename,mask_obj)

fnpattern_3dvols=d2n2s_write(mask_obj,tempdir,['3d_vol_write_' dicomuid],'vol',3);
mkdir(fileparts(outbasename));
outmat=[outbasename '.mat'];
outcon=[outbasename '.con'];


[~,outmsg]=system(['setup_masks ' inmat ' ' incon ' ' outbasename ' ' fnpattern_3dvols])

% get a randomise command addendum (the vxf and vxl inputs)
ind=strfind(outmsg,'--vxf=');
vxf_addendum=strrep(outmsg(ind:ind+7),' ','') %this will fail if mask dim 4 is 100 or more entries; this is a hitherto unencountered case so if you come across it think of something a little more clever for this section, big lol 

randomise_cmd_addendum=[' --vxf=' outbasename ' ' vxf_addendum];

end