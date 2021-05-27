
| Written by: | Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL                         |
| ----------- | ---------------------------------------------------------------------------------|
| When:       | February 2020 - May 2021                                                         |
| Purpose:    | Study "Beyond BOLD: in search of genuine diffusion fMRI contrast in human brain" |
| Paper:      | https://www.biorxiv.org/content/10.1101/2021.05.16.444253v1.full.pdf             |
| Contact:    | olszowyw@gmail.com                                                               |

Repeat the analyses
==============

We can share the datasets upon request.

Adjust the paths in the codes. Paths are only defined in the very first lines.

The order in which the scripts should be run from within MATLAB:
```
run('prepare_niftis.m');
run('make_VOI_masks.m');
run('fMRI_processing_single_subject.m');
run('GLM_single_subject_one_basis_function_only.m');
run('plot_individual_responses.m');
run('merge_figs_from_individuals.m');
run('summary_tSNRs.m');
run('group_analyses_task_spatial.m');
run('group_analyses_task_temporal.m');
run('group_analyses_rest.m');
```

Software
==============

I used the following software:

- Ubuntu 20.04.1
- MATLAB 2020a
- FSL 6.0.3
- SPM12 (r7771)
- ANTs 2.3.4

Repository contents
==============

- `paradigms`
   - `checkerboard_basic.py`
   
     PsychoPy code for the task.
   - `fixation_cross_3T.png`
   
     Fixation cross for resting-state, shown at the Prisma (3T).
   - `fixation_cross_7T.png`
   
     Fixation cross for resting-state, shown at the 7T.
- `software`

  Folder with packages and auxiliary functions used in the project.
   - `brain_parcellation`
   
      Files related to brain parcellation: atlases and masks.
   - `extra_functions`

      Auxiliary functions.
      - `plot_diagnostic_plots.m`
   
        Function to plot diagnostic plots: power spectra and QQ-plots.
      - `print_to_svg_to_pdf.m`

        Function to print a MATLAB figure to '.svg' and then to '.pdf', cropping the 'pdf' (removing margins) and deleting the '.svg' file at the end.
      - `spm_fMRI_design.m`
	  
		When loaded after SPM is loaded, this function overwrites the original SPM function: needed for VOI/PCA and FIR analyses.
      - `spm_regions.m`
	  
		When loaded after SPM is loaded, this function overwrites the original SPM function: needed for VOI/PCA and FIR analyses.
   - `mppca_denoise`
   
     Downloaded repository with an implementation of the MP-PCA denoising. For license, check `README.md` in that folder.
	 
   - `registration_to_MNI`
   
     Files needed for ANTs registrations.
      - `MICCAI2012-Multi-Atlas-Challenge-Data`
   
        Folder with files needed for ANTs brain extraction. For license, check https://figshare.com/articles/dataset/ANTs_ANTsR_Brain_Templates/915436
      - ...
   - `reisert-unring`
   
     Downloaded repository with an implementation of the Gibbs unringing. For license, check `/matlab/ringRm.cpp` in that folder.
   - `spm12`
   
     Downloaded SPM package. For license, check `LICENCE.txt` in that folder.
   - `acqp_AP.txt`
   
     Settings for topup if anterior-posterior phase encoding.
   - `acqp_LR.txt`
   
     Settings for topup if left-right phase encoding.
   - `design_rest.fsf`
   
     Settings for FSL MELODIC analysis of the resting-state data.
   - `design_task.fsf`
   
     Settings for FSL FEAT GLM analysis of the task data.
   - `mseq_stim_onsets.txt`

     Stimulus onset times for an m-sequence. Used in a pilot study.
   - `topup.cnf`
   
     Settings for topup. Taken from the UK BioBank MRI processing pipeline.
   
- `convert_pdfs_to_pngs.m`

  The other codes produce figures in pdf. They are saved in the 'figures' folder. This code produces pngs from these pdfs.
- `fMRI_processing_single_subject.m`

  Processing pipeline for SE BOLD and dfMRI datasets. Single-subject analyses only. Both task and resting state.
- `GLM_single_subject_one_basis_function_only.m`

  Running GLM analyses with only one basis function: the rectangular or canonical HRF.
- `group_analyses_rest.m`

  Performing resting-state (RS) group analyses: plotting group-averaged functional connectivity (FC) matrices, as well as running statistical tests.
- `group_analyses_task_spatial.m`

  Performing spatial task group analyses: averaging z-statistic maps across subjects. Looking how much activity there is in different ROIs.
- `group_analyses_task_temporal.m`

  Performing temporal task group analyses: plotting response functions. This analysis is run for voxels selected based on their z-statistics (transformed from F-statistics).
- `make_VOI_masks.m`

  Making VOI (= ROI) masks. Using NMM and JHU WM atlases. Needs to be run only once - before the single-subject analyses.
- `merge_figs_from_individuals.m`

  Merging smaller figures into bigger ones. Also, merging figures across subjects, producing long pdfs - with many pages, where each page refers to one subject/dataset.
- `omit_imgs_with_little_activation.m`

  Checking if enough activation in an ROI to consider the dataset in a group analysis. For task data only.
- `plot_individual_responses.m`

  Plotting mean time-courses and response functions for one exemplary subject at 3T and for one exemplary subject at 7T.
- `prepare_niftis.m`

  Preparing the NIfTI datasets for later analyses with `fMRI_processing_single_subject.m`.
- `summary_tSNRs.m`

  Summarising tSNR maps.
