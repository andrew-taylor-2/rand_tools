function make_config

config_fn=fullfile(fileparts(mfilename('fullpath')),'config.mat');

%% give prompts
full_functional_image_fn=input('give full filename for functional time series nifti');
full_data_coverage_image_fn=input('give full filename for data coverage/lesion time series nifti (to be used as a "lesion mask" in randomise)');
consensus_mask_fn=input('give full filename for consensus mask (where all time series has data -- not used if partial data option is true)');

%% check inputs -- file existence and try d2n2s on it
%functional
if ~exist(full_functional_image_fn,'file')
    error('couldn''t find the full functional time series file')
else
    try 
        d2n2s(full_functional_image_fn,'no','bvalbvecjson');
    catch 
        error('opening full functional time series file in matlab using spm failed')
    end
end

if ~exist(full_data_coverage_image_fn,'file')
    error('couldn''t find the full data coverage time series mask file')
else
    try 
        d2n2s(full_data_coverage_image_fn,'no','bvalbvecjson');
    catch 
        error('opening full data coverage time series mask file in matlab using spm failed')
    end
end

if ~exist(consensus_mask_fn,'file')
    error('couldn''t find the consensus mask file')
else
    try 
        d2n2s(consensus_mask_fn,'no','bvalbvecjson');
    catch 
        error('opening consensus mask file in matlab using spm failed')
    end
end

%% save
delete(config_fn)
save(config_fn)