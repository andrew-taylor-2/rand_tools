function [matoutfn,conoutfn,out1,out2]=make_design(mat_mat,con_mat,mat_fn,con_fn)
%write text file
dlmwrite(mat_fn,mat_mat,' ')

out1={};
out2={};
%turn to .mat
matoutfn=strrep(mat_fn,'.txt','.mat');
[out1{end+1},out2{end+1}]=system(['Text2Vest ' mat_fn ' ' matoutfn])

%write text file
dlmwrite(con_fn,con_mat,' ')

%turn to .con
conoutfn=strrep(con_fn,'.txt','.con');
[out1{end+1},out2{end+1}]=system(['Text2Vest ' con_fn ' ' conoutfn])
