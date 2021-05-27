

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Running GLM analyses with only one basis function: the rectangular or canonical HRF.
%%%%   Doing it, as the use of FIR basis functions might sometimes decrease sensitivity
%%%%   (more regressors). Run this code after performing single-subject analyses with
%%%%   'fMRI_processing_single_subject.m'.
%%%%
%%%%   I tested this code on one dataset (SE BOLD):
%%%%   BB_009_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200916093315_18
%%%%   with FIR 30/7 (no FAST). The results were similar to the FSL's output from the pipeline,
%%%%   though clearer for SPM than FSL: zfstat higher in visual cortex and lower elsewhere.
%%%%   Perhaps SPM has better microtime settings. Or due to pre-whitening differences.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage          = '/home/wiktor/Desktop/Dfmri/pipeline';
path_manage_software = fullfile(path_manage, 'software');


for analysis_id = 1:2
   
   if analysis_id == 1
      bf = 'rectangular';
   else
      bf = 'canonical';
   end
   
   folder_output = ['SPM_task_GLM_subject_space_one_' bf '_basis_function'];
   
   %-SPM
   addpath(genpath(fullfile(path_manage_software, 'spm12')));
   
   %-some auxiliary functions
   addpath(fullfile(path_manage_software, 'extra_functions'));
   
   %-dfMRI and SE BOLD images are 'adc' and 'ardfMRI_denoised_unringed_sdc', respectively
   task_folders = dir('/home/wiktor/Desktop/Dfmri/pipeline/analyses_output/*task*/a*.feat');
   
   for folder_id = 1:length(task_folders)
      
      cd(task_folders(folder_id).folder);
      
      if exist(fullfile(pwd, folder_output), 'dir') == 7
         continue
      end
      
      if strcmp(task_folders(folder_id).name, 'adc.feat')
         dataset_name = 'adc';
      else
         dataset_name = 'ardfMRI_denoised_unringed_sdc';
      end
      
      header_4D    = niftiinfo(dataset_name);
      No_volumes   = header_4D.ImageSize(4);
      TR           = header_4D.PixelDimensions(4);
      
      exper_design = readtable('exper_design.txt', 'Delimiter', ' ', 'ReadVariableNames', false);
      exper_design = table2array(exper_design);
      
      system(['gunzip ' dataset_name '.nii.gz']);
      
      scans = strings(No_volumes, 1);
      
      for volume_id = 1:No_volumes
         
         scans{volume_id} = [fullfile(pwd, [dataset_name '.nii']) ',' num2str(volume_id)];
         
      end
      
      clear jobs;
      
      %-model specification and estimation
      jobs{1}.stats{1}.fmri_spec.dir          = cellstr(fullfile(pwd, folder_output));
      jobs{1}.stats{1}.fmri_spec.timing.units = 'secs';
      jobs{1}.stats{1}.fmri_spec.timing.RT    = TR;
      
      %-with standard 16/8 problems with the contrast for the ADC F-test
      jobs{1}.stats{1}.fmri_spec.timing.fmri_t  = 32;
      jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 17;
      
      jobs{1}.stats{1}.fmri_spec.sess.scans     = cellstr(scans);
      
      %-specifying the experimental design
      jobs{1}.stats{1}.fmri_spec.sess.cond(1).name     = 'checkerboard';
      jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = exper_design(:,1);
      jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = exper_design(:,2);
      
      if strcmp(bf, 'rectangular')
         jobs{1}.stats{1}.fmri_spec.bases.none         = true;
      else
         jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs   = [0 0];
      end
      
      %-following https://www.nature.com/articles/s41467-019-09230-w.pdf
      jobs{1}.stats{1}.fmri_spec.cvi                   = 'FAST';
      
      %-no brain mask calculated
      jobs{1}.stats{1}.fmri_spec.mthresh   = -Inf;
      jobs{1}.stats{1}.fmri_spec.mask      = {''};
      
      jobs{1}.stats{2}.fmri_est.spmmat(1)  = cellstr([folder_output '/SPM.mat']);
      
      %-save residuals
      jobs{1}.stats{2}.fmri_est.write_residuals = 1;
      
      jobs{1}.stats{3}.con.spmmat          = cellstr([folder_output '/SPM.mat']);
      jobs{1}.stats{3}.con.consess{1}.tcon = struct('name', 'activation', 'convec', 1, 'sessrep', 'none');
      
      spm_jobman('run', jobs);
      
      %-merging SPM's GLM residuals to one file
      cd(folder_output);
      res_all    = dir('Res_*.nii');
      res_all    = {res_all.name};
      res_all_4d = char(res_all);
      
      %-TR can not be specified...
      spm_file_merge(res_all_4d, 'res4d.nii');
      for res_id = 1:No_volumes
         delete(char(res_all(res_id)));
      end
      
      cd ..
      
      plot_diagnostic_plots(fullfile(pwd, folder_output), 'res4d', 'none', TR, 1/128);
      
      cd(folder_output);
      
      load SPM;
      
      %-transforming the t-statistic map to a z-statistic map
      %-https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;4c24cff5.0808
      system('fslmaths spmT_0001.nii -mul 0 -add 1 ones.nii');
      system(['ttoz ones spmT_0001.nii ' num2str(SPM.xX.erdf) ' -zout zstat']);
      
      cd ..
      
      %-bringing zstat to MNI space
      
      image_name = pwd;
      image_name = image_name(strfind(image_name, 'BB_'):end);
      
      %-identifying the structural (MP2RAGE or MPRAGE) image
      anat_image    = dir(fullfile(path_manage, 'scans', [image_name(1:9) '*mp2rage*']));
      if isempty(anat_image)
         anat_image = dir(fullfile(path_manage, 'scans', [image_name(1:9) '*mprage*']));
      end
      anat_name = anat_image(1).name;
      anat_name = anat_name(1:end-7);
      
      anat2mni_affine_dir   = dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain0GenericAffine.mat']));
      anat2mni_warp_dir     = dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain1Warp.nii.gz']));
      anat2mni_inv_warp_dir = dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain1InverseWarp.nii.gz']));
      anat2mni_affine       = fullfile(anat2mni_affine_dir(1).folder,   anat2mni_affine_dir(1).name);
      anat2mni_warp         = fullfile(anat2mni_warp_dir(1).folder,     anat2mni_warp_dir(1).name);
      anat2mni_inv_warp     = fullfile(anat2mni_inv_warp_dir(1).folder, anat2mni_inv_warp_dir(1).name);
      
      system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
         ' -i ' fullfile(pwd, folder_output, 'zstat.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
         ' -o ' fullfile(pwd, folder_output, 'zstat_MNI.nii.gz')]);
      
      system(['gzip ' dataset_name '.nii']);
      
   end
   
end
