%% Sample script to call eddy current correction functions in FSL

% FSL root directory
FSL_dir = '/usr/local/fsl';

% Prepare parameter structure
study_info.label1 = 'Subj09';
study_info.label2 = 'active';
study_info.label3 = 's1';
study_info.label4 = 'NSA4';
study_info.PE = [0, 1, 0];
study_info.PE_aux = [0, -1, 0];
study_info.b0_idx = 24;
study_info.b0_idx_aux = 0;
study_info.readout = 0.8;
study_info.readout_aux = 0.8;
study_info.mask_thres = 0.6;
study_info.mask_tag = 'muscle';
study_info.scan_dir = 'my/input_dir';
study_info.dwi_tag = 'my_DWI_img';
study_info.dwi_tag_aux = 'my_DWI_img_b0';
study_info.output_dir = 'my/output_dir';

% Perform eddy current correction
eddy_correct(FSL_dir, study_info);
