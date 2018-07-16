This folder contains custom made scripts and example dataset which demonstrate SPM longitudinal processing 
of qMRI data applied during this study. 

1. Install software
Before beeing able to use this code please install SPM and include in matlab path (https://www.fil.ion.ucl.ac.uk/spm/software/download/)
Please also install CAT toolbox for SPM into toolbox subdirectory 'cat12' (http://www.neuro.uni-jena.de/cat/index.html#DOWNLOAD);
Finally, please copy the "trajectory" folder into the SPM toolbox subdirectory before continue.

2. unzip image files in all “mpm_example_data” subfolders

3. To run Voxel-based Morphometry (VBM) preprocessing demo with an example dataset of one anonymous subject
please run matlab script "demo_preprocessing_vbm"

4. To run Voxel-based Quantification (VBQ) preprocessing demo with all MPM parameter maps of the same subject
please run "demo_preprocessing_vbq"

The expected runtime based on a MACBOOK PRO 2,5 GHz Intel Core i7 is around an 1-2 hours. 

Please note that this particular code has not been thoroughly tested with other software versions than MATLAB2016b, SPM12 r7355, and CAT12 r1318 and
other parameter choices might produce errors during processing attempts. The code aims at transparency and illustration but is not intended for clinical use.

The ./covariate subfolder contains a so called person period variable in (PP.mat) with information about the sample 
(here only 1 subject with 2 scans but his can be easily extended).

The ./mpm_example_data/mpm_generated_nativespace subfolder contains MPM raw data of one subject with two scans including MT, R1, R2* and PD parameters 
(obtained from MPM map generation module using hMRI toolbox for SPM).

The ./mpm_example_data/mpm_spatially_processed subfolder contains previously calculated processed MPM data of the example subject.

The ./trajectory subfolder contains a reduced version of custom SPM toolbox code with two functions performing VBM and VBQ respectively. 
After installing the toolbox in SPM (ie. copying it to SPM toolbox folder) it can be used in the usual way with the batch job manager (spm_jobman)
or the SPM GUI.

This is free but copyright software, distributed under the terms of the GNU General Public Licence as published by the Free Software Foundation (either version 2, or at your option, any later version). Further details on "copyleft" can be found at http://www.gnu.org/copyleft/. In particular, software is supplied as is. No formal support or maintenance is provided or implied.
  
For any questions and requests please contact gabriel.ziegler@dzne.de

