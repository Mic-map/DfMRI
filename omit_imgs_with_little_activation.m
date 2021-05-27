

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Checking if enough activation in an ROI to consider the dataset in a group analysis. Run
%%%%   from within 'make_group_analyses_task_temporal.m'.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       August 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd(fullfile(path_manage, 'analyses_output', imgs(img_id)));

omit       = false;

mask_VOI   = niftiread(['mask_VOI_' VOI '_subject_space_dilated.nii.gz']);

if     subplot_id == 1 || subplot_id == 6
   zfstat1 = niftiread(fullfile('ardfMRI_denoised_unringed_sdc.feat', 'stats', 'zfstat1.nii.gz'));
elseif subplot_id == 3 || subplot_id == 8
   zfstat1 = niftiread(fullfile('map_0pt2.feat',                      'stats', 'zfstat1.nii.gz'));
else
   zfstat1 = niftiread(fullfile('adc.feat',                           'stats', 'zfstat1.nii.gz'));
end

%-arbitrary threshold
mask_VOI(zfstat1 < 3.1) = 0;

if sum(reshape(mask_VOI, 1, [])) < 20
   omit    = true;
end
