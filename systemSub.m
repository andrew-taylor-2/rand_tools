
function status  = systemSub (cmd)
% Save library paths
MatlabPath = getenv('LD_LIBRARY_PATH');
% Make Matlab use system libraries
setenv('LD_LIBRARY_PATH',getenv('PATH'));
fprintf('%s\n',cmd);
status = system(cmd);
% Reassign old library paths
setenv('LD_LIBRARY_PATH',MatlabPath);
%systemSub