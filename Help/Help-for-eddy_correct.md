# Help for the function <i>eddy_correct</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#1-purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#2-usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#4-Example-Code)
5) [Useful Links](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#5-Useful-Links)
6) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#6-Acknowledgements)


## 1. Purpose

The function <i>eddy_correct</i> is used to perform the eddy current correction on DWI images by invoking FSL built-in tools [<i>topup</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup) and [<i>eddy</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy).

## 2. Usage
This function must be run under Unix-based systems. To perform eddy current correction, the user should have all images to be processed in NIfTI format.

The user provides the root directory of FSL and a structure containing all required parameters (see [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#3-Syntax) for detailed description). The function first corrects the susceptibility induced distortions by invoking [<i>topup</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup) from FSL, and then performs eddy current correction by invoking [<i>eddy</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy).

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md)

## 3. Syntax
eddy_correct(FSL_dir, DTI_scan_info)

The input arguments are:
 * FSL_dir: A character string specifying the root directory for FSL
 * DTI_scan_info: A structure containing information of the images
   * subject: A character string specifying the subject (e.g. 'Subj09')
   * status: A character string specifying the subject conditions under which the images are acquired (e.g. 'active', 'passive', etc.)
   * stack: A character string specifying the stack label of the scan
   * scan_label: A character string for additional label
   * PE: A 1-by-3 vector specifying the phase-encoding direction [x, y, z]; some commonly applied options are:
 
     P->A: [0,  1, 0]

     A->P: [0, -1, 0]

     L->R: [-1, 0, 0]

     R->L: [1,  0, 0]
 
   * PE_aux: Phase-encoding direction for auxiliary b=0 image
   * b0_idx: Frame index with b=0 in DWI image. **Note that FSL treats index starting from 0.**
   * b0_idx_aux: Frame index with b=0 in auxiliary image
   * readout: Total readout time (sec). Refer to [<i>FSL wiki</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide#A--datain) for more details.
   * readout_aux: Total readout time for auxiliary b=0 image
   * mask_thres: Fraction intensity threshold for binary mask
   * mask_tag: A character string specifying the name tag for the binary mask (e.g. 'muscle')
   * scan_dir: A character string specifying the directory where raw image data are stored
   * dwi_tag: A character string specifying the main DWI image file name **without the '.nii' extention**
   * dwi_tag_aux: A character string specifying the auxiliary (b=0) DWI image file name **without the '.nii' extention**
   * output_dir: A character string specifying the root directory to store all processed images and data

There is no return value from this function. For the directories where the images are read and processed data saved, refer to [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md#4-Example-Code) below.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md)

## 4. Example Code
    % FSL root directory
    FSL_dir = '/usr/local/fsl';

    % Prepare parameter structure
    study_info.subject     = 'Subj09';
    study_info.status      = 'active';
    study_info.stack       = 's1';
    study_info.scan_label  = 'NSA4';
    study_info.PE          = [0, 1, 0];
    study_info.PE_aux      = [0, -1, 0];
    study_info.b0_idx      = 24;
    study_info.b0_idx_aux  = 0;
    study_info.readout     = 0.8;
    study_info.readout_aux = 0.8;
    study_info.mask_thres  = 0.6;
    study_info.mask_tag    = 'muscle';
    study_info.scan_dir    = 'my/input_dir';
    study_info.dwi_tag     = 'my_DWI_img';
    study_info.dwi_tag_aux = 'my_DWI_img_b0';
    study_info.output_dir  = 'my/output_dir';

    % Perform eddy current correction
    eddy_correct(FSL_dir, study_info);

In this example, the function reads image files <i>my/input_dir/my_DWI_img.nii</i> and <i>my/input_dir/my_DWI_img_b0.nii</i> and keeps a copy of unprocessed images in directory <i>my/output_dir/Subj09/active/s1/raw</i>. All processed data by <i>topup</i> and <i>eddy</i> are stored in directory <i>my/output_dir/Subj09/active/s1/eddy</i>

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md)

## 5. Useful Links
 * [A general tutorial for FSL Diffusion Toolbox](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/fdt1/index.html)
 * [FSL wiki for <i>topup</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup)
 * [FSL wiki for <i>eddy</i>](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy)

## 6. Acknowledgements
Grants:

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-eddy_correct.md)
