function eddy_correct(FSL_dir, DTI_scan_info)
%
% FUNCTION eddy_correct
%
% USAGE
%   Perform eddy current correction through FSL
%
% INPUT ARGUMENTS
%   FSL_dir          : FSL directory
%   DTI_scan_info    : Structure specifying DTI image information
%       -label1      : First label for the image file names, e.g. 'S09'
%       -label2      : Second label for the image file names, e.g. 'active' or 'passive'
%       -label3      : Third label for the image file names, e.g. 's1' or 's2'; lower case
%       -label4      : Addition label for the scan, e.g. 'NSA4_10'
%       -PE          : Phase-encoding direction; [x, y, z]
%       -PE_aux      : Phase-encoding direction for auxiliary b=0 image
%       -b0_idx      : Frame index with b=0 in DWI image
%       -b0_idx_aux  : Frame index with b=0 in auxiliary image
%       -readout     : Total readout time (sec)
%       -readout_aux : Total readout time for auxiliary b=0 image
%       -mask_thres  : Fraction intensity threshold for binary mask
%       -mask_tag    : Tag string for binary mask, e.g. 'muscle'
%       -scan_dir    : Directory where raw image data are stored
%       -dwi_tag     : Main DWI image file name prefix (w.o. '.nii')
%       -dwi_tag_aux : Auxiliary (b=0) DWI image file name prefix
%       -output_dir  : Root directory for all processed images
%
% VERSION INFOMRATION
%   v. 1.0.0 Xingyu Zhou 

%% 00 - Initialization
FSL_env_setup(FSL_dir);

% Directory for processed subject data
subj_out_dir = sprintf('%s/%s/%s/%s', ...
    DTI_scan_info.output_dir, ...
    DTI_scan_info.label1, ...
    DTI_scan_info.label2, ...
    DTI_scan_info.label3);

if ~exist(subj_out_dir, 'dir')
    mkdir(subj_out_dir);
end

%% 01.0 - Prepare directories and image files to read
raw_dir = sprintf('%s/%s', subj_out_dir, 'raw');
eddy_dir = sprintf('%s/%s', subj_out_dir, 'eddy');

if ~exist(raw_dir, 'dir')
    mkdir(raw_dir);
end

if ~exist(eddy_dir, 'dir')
    mkdir(eddy_dir);
end

fdwi    = sprintf('%s/%s', ...
    DTI_scan_info.scan_dir, DTI_scan_info.dwi_tag);
fb0_aux = sprintf('%s/%s', ...
    DTI_scan_info.scan_dir, DTI_scan_info.dwi_tag_aux);
img_tag = sprintf('%s_%s_%s_%s', ...
    DTI_scan_info.label1, ...
    DTI_scan_info.label2, ...
    DTI_scan_info.label3, ...
    DTI_scan_info.label4);
raw_base_str = sprintf('%s/%s', raw_dir, img_tag);
eddy_base_str = sprintf('%s/%s', eddy_dir, img_tag);

dwi_info = niftiinfo([fdwi '.nii']);
img_size = dwi_info.ImageSize;
clear dwi_info;

%% TOPUP 00 - Trim images into even x, y, z dimensions
% Display basic image info using FSLINFO
fprintf('\n=== Basic image info of DWI image ===\n  %s\n', fdwi);
system(['fslinfo ' fdwi]);
fprintf('\n=== Basic image info of auxiliary b = 0 image ===\n  %s\n\n', fb0_aux);
system(['fslinfo ' fb0_aux]);

fdwi_raw = sprintf('%s_%s', raw_base_str, 'DWI');
fb0_aux_raw = sprintf('%s_%s', raw_base_str, 'b0_aux');
cmd_trim_dwi = strjoin({'fslroi', ...
                        ['"' fdwi '"'], ...
                        ['"' fdwi_raw '"'], ...
                        '0', int2str(img_size(1) - mod(img_size(1), 2)), ...
                        '0', int2str(img_size(2) - mod(img_size(2), 2)), ...
                        '0', int2str(img_size(3) - mod(img_size(3), 2))});
cmd_trim_b0_aux = strjoin({'fslroi', ...
                        ['"' fb0_aux '"'], ...
                        ['"' fb0_aux_raw '"'], ...
                        '0', int2str(img_size(1) - mod(img_size(1), 2)), ...
                        '0', int2str(img_size(2) - mod(img_size(2), 2)), ...
                        '0', int2str(img_size(3) - mod(img_size(3), 2))});
status = system(cmd_trim_dwi);
if status == 0
    fprintf('DWI raw image saved to even x, y, z dimensions.\n');
end

status = system(cmd_trim_b0_aux);
if status == 0
    fprintf(['Auxiliary b = 0 raw image saved to ' ...
             'even x, y, z dimensions.\n']);
end

%% TOPUP 01.1 - Extract b = 0 frame from DWI image set using fslroi
fb0_main_save = sprintf('%s_%s', eddy_base_str, 'b0');
cmd_extract_b0 = strjoin({'fslroi', ...
                          ['"' fdwi_raw '"'], ...
                          ['"' fb0_main_save '"'], ...
                          int2str(DTI_scan_info.b0_idx), '1'});
status = system(cmd_extract_b0);
if status == 0
    fprintf('FSLROI : extracted b = 0 frame from raw DWI image.\n');
end

%% TOPUP 01.2 - Extract b = 0 frame from auxiliary image using fslroi
fb0_aux_save = sprintf('%s_%s', eddy_base_str, 'b0_aux');
cmd_extract_b0_aux = strjoin({'fslroi', ...
                          ['"' fb0_aux_raw '"'], ...
                          ['"' fb0_aux_save '"'], ...
                          int2str(DTI_scan_info.b0_idx_aux), '1'});
status = system(cmd_extract_b0_aux);
if status == 0
    fprintf(['FSLROI : extracted b = 0 frame ' ...
             'from raw auxiliary image.\n']);
end

%% TOPUP 02 - Merge b = 0 images into a single image along the 4th ("time") axis
fb0_merge = sprintf('%s_%s', eddy_base_str, 'b0_merge');
cmd_merge = strjoin({'fslmerge -t', ...
                     ['"' fb0_merge '"'], ...
                     ['"' fb0_main_save '"'], ...
                     ['"' fb0_aux_save '"']});
status = system(cmd_merge);
if status == 0
    fprintf(['FSLMERGE : b = 0 frames with different PE directions ' ...
             'merged along the 4th axis.\n']);
end

%% TOPUP 03 - Prepare parameter file
fparam = sprintf('%s_%s', eddy_base_str, 'acqparams.txt');
param_mat = [DTI_scan_info.PE, DTI_scan_info.readout; ...
             DTI_scan_info.PE_aux, DTI_scan_info.readout_aux];
dlmwrite(fparam, param_mat, ' ');
fprintf('Acquisition parameter file for TOPUP ready.\n');

%% TOPUP 04 - Prepare and execute ultimate command for topup
ftopup_base = sprintf('%s/%s_%s', eddy_dir, 'topup', img_tag);
ftopup_iout = sprintf('%s_%s', ftopup_base, 'iout');
ftopup_fout = sprintf('%s_%s', ftopup_base, 'fout');
fconfig     = 'b02b0.cnf';    % FSL built-in configuration file for topup

cmd_topup = strjoin({'topup', ['--imain="'  fb0_merge '"'], ...
                              ['--datain="' fparam '"'], ...
                              ['--config="' fconfig '"'], ...
                              ['--out="'    ftopup_base '"'], ...
                              ['--iout="'   ftopup_iout '"'], ...
                              ['--fout="'   ftopup_fout '"']});
fprintf('Running TOPUP...\n');
status = system(cmd_topup);
if status == 0
    fprintf('TOPUP : command finished successfully.\n');
end

%% EDDY 01 - Compute average image of corrected b = 0 volumes
favg = sprintf('%s_%s', eddy_base_str, 'b0_avg_cor');
cmd_avg = strjoin({'fslmaths', ...
                   ['"' ftopup_iout '"'], ...
                   '-Tmean', ...        % avg along "time" axis
                   ['"' favg '"']});
status = system(cmd_avg);
if status == 0
    fprintf(['FSLMATHS : computed average image of ' ...
             'corrected b = 0 volumes.\n']);
end

%% EDDY 02 - Create binary mask
fthres = sprintf('%s_%s_%s', eddy_base_str, 'b0', DTI_scan_info.mask_tag);
cmd_mask = strjoin({'bet', ...
                    ['"' favg '"'], ...
                    ['"' fthres '"'], ...
                    '-m', '-f', num2str(DTI_scan_info.mask_thres)});
status = system(cmd_mask);
if status == 0
    fprintf(['BET : binary mask generated; intensity threshold = ' ...
            num2str(DTI_scan_info.mask_thres) '\n']);
end

%% EDDY 03 - Prepare index file
findex = sprintf('%s_%s', eddy_base_str, 'index.txt');
dlmwrite(findex, ones(img_size(4), 1));
fprintf('Index file for EDDY ready.\n');

%% EDDY 04 - Prepare b-vector and b-value files
fbvecs_raw = sprintf('%s/%s%s', ...
    DTI_scan_info.scan_dir, DTI_scan_info.dwi_tag, '.bvec');
fbvals_raw = sprintf('%s/%s%s', ...
    DTI_scan_info.scan_dir, DTI_scan_info.dwi_tag, '.bval');

fbvecs = sprintf('%s%s', raw_base_str, '.bvec');
fbvals = sprintf('%s%s', raw_base_str, '.bval');

copyfile(fbvecs_raw, fbvecs);
fprintf('File for b-vectors ready.\n');
copyfile(fbvals_raw, fbvals);
fprintf('File for b-values ready.\n');

%% EDDY 05 - Prepare and execute ultimate command for eddy
fmask      = sprintf('%s_%s', fthres, 'mask');
fwhm       = 0;
flm        = 'quadratic';
feddy_base = sprintf('%s/%s_%s', eddy_dir, 'eddy_unwarped', img_tag);

cmd_eddy = strjoin({'eddy', ['--imain="' fdwi_raw '"'], ...
                            ['--mask="'  fmask '"'], ...
                            ['--index="' findex '"'], ...
                            ['--acqp="'  fparam '"'], ...
                            ['--bvecs="' fbvecs '"'], ...
                            ['--bvals="' fbvals '"'], ...
                            ['--fwhm="'  num2str(fwhm) '"'], ...
                            ['--topup="' ftopup_base '"'], ...
                            ['--flm="'   flm '"'], ...
                            ['--out="'   feddy_base '"'], ...
                            '--data_is_shelled'});
fprintf('Running EDDY...\n');
status = system(cmd_eddy);
if status == 0
    fprintf('EDDY : command finished successfully.\n');
end

end

function FSL_env_setup(FSL_dir)
%FSL_ENV_SETUP Initialize FSL environment variables
%   $FSLDIR        : FSL root directory
%   $PATH          : System PATH
%   $FSLOUTPUTTYPE : FSL output image format (.nii.gz)

FSL_path = sprintf('%s/%s', FSL_dir, 'bin');

if ~contains(getenv('FSLDIR'), FSL_dir)
    setenv('FSLDIR', FSL_dir);
    fprintf('\nEnvironment variable FSLDIR set.\n');
end

if ~contains(getenv('PATH'), FSL_path)
    setenv('PATH', [getenv('PATH') ':' FSL_path]);
    fprintf('\nEnvironment variable PATH set.\n');
end

if ~contains(getenv('FSLOUTPUTTYPE'), 'NIFTI_GZ')
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    fprintf('\nEnvironment variable FSLOUTPUTTYPE set.\n');
end

end
  
