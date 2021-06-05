# README.md
This file has information about the sample data files.

## File named data_for_fibertracking.mat
This file has the selected outputs from [pre_process.m](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Sample-Scripts/pre_process.m).  When the file is loaded, type

  whos

on the command line. The output is:

  Name   Size   Bytes   Class   Attributes

  anat_fov  1x2   16  double              
  
  anat_image           192x192x44             12976128  double              
  
  anat_numcols           1x1                         8  double              
  
  anat_numrows           1x1                         8  double              
  
  anat_numslcs           1x1                         8  double              
  
  anat_pixelspace        1x2                        16  double              
  
  anat_slcthick          1x1                         8  double              
  
  dti_all_numslcs        1x1                         8  double              
  
  dti_fov                1x2                        16  double              
  
  dti_numcols            1x1                         8  double              
  
  dti_numdir             1x1                         8  double              
  
  dti_numrows            1x1                         8  double              
  
  dti_pixelspace         1x2                        16  double              
  
  dti_slcthick           1x1                         8  double              
  
  tensor_m               5-D                 116785152  double              

These variables are needed to perform the fiber-tract processing included in [process_tracts.m](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Sample-Scripts/process_tracts.m).

## File named ta_mask.mat

This file has the muscle mask used in the manuscript.  After loading it, type

whos

The output is: 

  Name        Size                   Bytes  Class     Attributes

  mask      192x192x44            12976128  double  
  
  
## File named ta_mesh.mat

This file has the aponeurosis mesh used in the manuscript.  After loading it, type

whos

  Name                      Size                   Bytes  Class     Attributes

  roi_mask                192x192x44            12976128  double              
 
  roi_mesh                182x36x6                314496  double              
  
  roi_mesh_dilated_1      182x36x6                314496  double              
