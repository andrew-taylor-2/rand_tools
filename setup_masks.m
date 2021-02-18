function [outmat,outcon,outmsg]=setup_masks(inmat,incon,outbasename,mask_obj)

fnpattern_3dvols=d2n2s_write(mask_obj,tempdir,['3d_vol_write_' dicomuid],'vol',3);
mkdir(fileparts(outbasename));
outmat=[outbasename '.mat'];
outcon=[outbasename '.con'];


[~,outmsg]=system(['setup_masks ' inmat ' ' incon ' ' outbasename ' ' fnpattern_3dvols])
end