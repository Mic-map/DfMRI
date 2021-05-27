

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Processing pipeline for SE BOLD and dfMRI datasets. Single-subject analyses only.
%%%%   Both task and resting state. In the variable names, 'b_0' denotes the lower b-value,
%%%%   while 'b_1' denotes the higher b-value. During the pilot study, we acquired data for
%%%%   b-values 0/1, thus, the naming. For a subject in one session (3T or 7T), it is important
%%%%   to run this pipeline in order (across 'analysis_id'): for a 'revPE' dataset, 'topup' runs
%%%%   only with the best aligned 'linPE' image, for the other 'linPE' datasets, 'applytopup'
%%%%   runs with the output of the saved 'topup'. If the pipeline runs on the datasets of the
%%%%   subject in a different order than following 'analysis_id', then susceptibility distortion
%%%%   correction (SDC) can be suboptimal.
%%%%
%%%%   It is assumed that all relevant datasets are in '/pipeline/scans/', and that the necessary
%%%%   initial preprocessing was done with 'prepare_niftis.m'. If the PE alignment is different
%%%%   than AP for RS and LR for task, then add an exception under '%-topup exceptions'.
%%%%
%%%%   It is crucial that 'No_dummies_out' is well defined - for that, always look at the NIfTIs
%%%%   before running this code - the first volumes with b-values 0 are usually redundant, and
%%%%   this needs to be later accounted for in the processing. If there are 3 volumes in the
%%%%   dataset which are redundant for later analyses, define 'No_dummies_out' as 3. Put it just
%%%%   after you define 'image_name' and 'revPE' - names of the images (without the extension).
%%%%
%%%%   Due to RAM problems, I remove many variables once they are not needed - with 'clear'.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       February 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LASTN                  = maxNumCompThreads(10);
path_manage            = '/home/wiktor/Desktop/Dfmri/pipeline';
path_manage_software   = fullfile(path_manage, 'software');
replace_dim_estimation = '0';
replace_output_comps   = '40';
frac_outliers_discard  = 0.2;


%-SPM
addpath(genpath(fullfile(path_manage_software, 'spm12')));

%-PCA denoising:     codes of Jelle Veraart
addpath(genpath(fullfile(path_manage_software, 'mppca_denoise')));

%-Gibbs unringing:   codes from Freiburg
addpath(genpath(fullfile(path_manage_software, 'reisert-unring-ddd43da65219')));

%-some auxiliary functions
%-needs to be added after SPM! as I changed 'spm_regions.m'
addpath(fullfile(path_manage_software, 'extra_functions'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-reading the Neuromorphometrics atlas parcellation

%-when 'readtable' is used, running later 'unring' leads to a mex file error
%-this is very surprising!
%-'textscan' is a way around
fileID = fopen(fullfile(path_manage_software, 'brain_parcellation', 'labels_Neuromorphometrics.txt'));
A      = textscan(fileID, '%s', 'Delimiter', '', 'headerLines', 1);
A      = A{1};

ROIs   = zeros  (1, size(A, 1));
ROIs_t = strings(1, size(A, 1));
for ROI_id  = 1:size(A, 1)
   full_ROI = A{ROI_id};
   ROIs(ROI_id)   = str2double(extractBetween(full_ROI, '<index>', '</index'));
   ROIs_t(ROI_id) = extractBetween(full_ROI, '<name>',  '</name>');
   %-shortening ROI names
   ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), 'Left',   'L');
   ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), 'Right',  'R');
   ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), 'gyrus',  'G');
   ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), 'cortex', 'C');
   ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), ' the',   '');
   ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), 'cortex', 'C');
   %-removing abbreviations in ROI names
   for abbr = {'PT', 'FO', 'AnG', 'CO', 'OpIFG', 'SOG', 'POrG', 'PO', 'TTG', 'MOG', 'STG', 'SPL', 'PCu', 'SMG', 'TrIFG', 'PP', 'SMC', 'Cun', 'MTG', 'MFG', 'OCP', 'AIns', 'MSFG', 'SFG', 'IOG', 'MPoG', 'Calc', 'PIns', 'TMP', 'SCA', 'LOrG', 'PoG', 'MPrG', 'PrG', 'Ent', 'OrIFG', 'PCgG', 'MCgG', 'MFC', 'AOrG', 'ACgG', 'FRP', 'OFuG', 'PHG', 'ITG', 'FuG', 'LiG', 'GRe', 'MOrG'}
      ROIs_t(ROI_id) = strrep(ROIs_t(ROI_id), [' ' abbr{1}], '');
   end
end
ROIs_t           = cellstr(ROIs_t);

atlas_path       = fullfile(path_manage_software, 'brain_parcellation', 'rlabels_Neuromorphometrics.nii');
atlas_w_JHU_path = fullfile(path_manage_software, 'brain_parcellation', 'NMM_JHU_atlases.nii.gz');


for analysis_id  = 1:196
   
   %-default settings
   revPE            = '';
   b_val_diff       = 0;
   vol_outliers_thr = 0.01;
   
   %-as we wanted to cover different b-values (e.g. calibration issues), for some datasets, the very first volumes were for b-value 1000 or 0, so sometimes need to be removed as dummies
   No_dummies_out = 0;
   
   if analysis_id == 1
      image_name  = 'pilot_RS_cmrr_mbep2d_se_b0_task_12s_off_18s_on_MB1_20200217155721_9';  %-poor func2mni, but this dataset is not important (pilot); SDC not great
      revPE       = 'pilot_RS_cmrr_mbep2d_b0_revPE_20200217155721_8';
   elseif analysis_id == 2
      image_name  = 'pilot_RS_cmrr_mbep2d_diff_b0_1_task_12s_off_18s_on_MB1_20200217155721_7';  %-poor func2mni, but this dataset is not important (pilot); SDC not great
      revPE       = 'pilot_RS_cmrr_mbep2d_b0_revPE_20200217155721_8';
      %%%
   elseif analysis_id == 3
      image_name  = 'pilot_WO_Dfmri_64cx_2iso_diff_600dir_b0_1_MB2_rs_20200227182627_10';
   elseif analysis_id == 4
      image_name  = 'pilot_WO_Dfmri_64cx_2iso_se_600dir_b0_MB2_rs_20200227182627_11';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 5
      image_name  = 'pilot_EF_Dfmri_64cx_sms2_2iso_diff_600dir_b0_1_task_12s_off_18s_on_20200312162637_9';  %-there was poor registration to MNI space when using the T1 map as the brain mask was poor, so I chose a different anat map
   elseif analysis_id == 6
      image_name  = 'pilot_EF_Dfmri_64cx_sms2_2iso_se_600dir_task_12s_off_18s_on_20200312162637_8';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 7
      image_name  = 'pilot_WO2_cmrr_mbep2d_se_task_18s_off_12s_on_20200520100309_12';
      revPE       = 'pilot_WO2_cmrr_mbep2d_se_revPE_20200520100309_9';
   elseif analysis_id == 8
      image_name  = 'pilot_WO2_cmrr_mbep2d_diff_b0_1_task_18s_off_2s_on_20200520100309_7';
      revPE       = 'pilot_WO2_cmrr_mbep2d_se_revPE_20200520100309_9';
   elseif analysis_id == 9
      image_name  = 'pilot_WO2_cmrr_mbep2d_diff_b0_1_task_18s_off_12s_on_20200520100309_8';
      revPE       = 'pilot_WO2_cmrr_mbep2d_se_revPE_20200520100309_9';
      %%%
   elseif analysis_id == 10
      image_name  = 'pilot_SM_cmrr_mbep2d_se_task_18s_off_2s_on_20200520110644_7';
      %%%
   elseif analysis_id == 11
      image_name  = 'BB_001_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200611094346_15';
      revPE       = 'BB_001_7T_cmrr_mbep2d_se_revPE_20200611094346_17';
   elseif analysis_id == 12
      image_name  = 'BB_001_7T_cmrr_mbep2d_se_rs_20200611094346_27';
      revPE       = 'BB_001_7T_cmrr_mbep2d_se_rs_revPE_20200611094346_25';
   elseif analysis_id == 13
      image_name  = 'BB_001_7T_cmrr_mbep2d_diff_b0_1_bipolar_task_18s_off_12s_on_20200611094346_11';
      revPE       = 'BB_001_7T_cmrr_mbep2d_diff_b0_1_bipolar_revPE_20200611094346_13';
   elseif analysis_id == 14
      image_name  = 'BB_001_7T_cmrr_mbep2d_diff_b0_1_bipolar_rs_20200611094346_21';
      revPE       = 'BB_001_7T_cmrr_mbep2d_diff_b0_1_bipolar_rs_revPE_20200611094346_23';
      %%%
   elseif analysis_id == 15
      image_name  = 'BB_002_7T_cmrr_mbep2d_diff_b0_1_task_18s_off_12s_on_20200630150831_22';
      revPE       = 'BB_002_7T_cmrr_mbep2d_diff_b0_1_revPE_20200630150831_14';
   elseif analysis_id == 16
      image_name  = 'BB_002_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200630150831_24';  %-wrong alignment 'BB_002_7T_cmrr_mbep2d_se_revPE_20200630150831_16'
      revPE       = 'BB_002_7T_cmrr_mbep2d_diff_b0_1_revPE_20200630150831_14';  %-different number of slices than originally linPE!
   elseif analysis_id == 17
      image_name  = 'BB_002_7T_cmrr_mbep2d_diff_b0_1_rs_20200630150831_26';
      revPE       = 'BB_002_7T_cmrr_mbep2d_diff_b0_1_rs_revPE_20200630150831_18';
   elseif analysis_id == 18
      image_name  = 'BB_002_7T_cmrr_mbep2d_se_rs_20200630150831_28';
      revPE       = 'BB_002_7T_cmrr_mbep2d_diff_b0_1_rs_revPE_20200630150831_18';  %-different number of slices than originally linPE!
      %%%
   elseif analysis_id == 19
      image_name  = 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200702093313_18';  %-subject was not lying straight in the scanner?
      revPE       = 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_true_1_revPE_20200702093313_12';
   elseif analysis_id == 20
      image_name  = 'BB_003_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200702093313_20';
      revPE       = 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_true_1_revPE_20200702093313_12';  %-slightly different spaces! imperfect topup; no true revPE... 'BB_003_7T_cmrr_mbep2d_se_revPE_20200702093313_10'
   elseif analysis_id == 21
      image_name  = 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20200702093313_22';
      revPE       = 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200702093313_14';
   elseif analysis_id == 22
      image_name  = 'BB_003_7T_cmrr_mbep2d_se_rs_20200702093313_24';
      revPE       = 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200702093313_14';  %-wrong alignment 'BB_003_7T_cmrr_mbep2d_se_revPE_rs_20200702093313_16'
      %%%
   elseif analysis_id == 23
      image_name  = 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200706151255_16';
      revPE       = 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_20200706151255_8';
   elseif analysis_id == 24
      image_name  = 'BB_004_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200706151255_18';
      revPE       = 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_20200706151255_8';  %-wrong alignment 'BB_004_7T_cmrr_mbep2d_se_revPE_20200706151255_10'
   elseif analysis_id == 25
      image_name  = 'BB_004_7T_cmrr_mbep2d_diff_b0pt29_1pt45_task_18s_off_12s_on_20200706151255_20';
      revPE       = 'BB_004_7T_cmrr_mbep2d_diff_b0pt29_1pt45_revPE_20200706151255_12';  %-maybe wrong alignment
      %%%
   elseif analysis_id == 26
      image_name  = 'BB_003_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200710180104_8';
      revPE       = 'BB_003_3T_cmrr_mbep2d_diff_0pt2_1_revPE_20200710180104_9';
   elseif analysis_id == 27
      image_name  = 'BB_003_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200710180104_11';
      revPE       = 'BB_003_3T_cmrr_mbep2d_se_revPE_20200710180104_15';
   elseif analysis_id == 28
      image_name  = 'BB_003_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200710180104_13';
      revPE       = 'BB_003_3T_cmrr_mbep2d_se_revPE_20200710180104_15';
   elseif analysis_id == 29
      image_name  = 'BB_003_3T_cmrr_mbep2d_diff_0pt29_1pt45_task_18s_off_12s_on_20200710180104_16';
      revPE       = 'BB_003_3T_cmrr_mbep2d_diff_0pt29_1pt45_revPE_20200710180104_17';
   elseif analysis_id == 30
      image_name  = 'BB_003_3T_cmrr_mbep2d_diff_0pt2_1_TR350_task_18s_off_12s_on_20200710180104_18';
      revPE       = 'BB_003_3T_cmrr_mbep2d_diff_0pt2_1_TR350_revPE_20200710180104_19';
      %%%
   elseif analysis_id == 31
      image_name  = 'BB_005_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200731100958_14';
      revPE       = 'BB_005_7T_cmrr_mbep2d_se_revPE_20200731100958_10';
   elseif analysis_id == 32
      image_name  = 'BB_005_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200731112951_5';
      revPE       = 'BB_005_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_20200731112951_3';
      %%%
   elseif analysis_id == 33
      image_name  = 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_rs_20200804162733_9';
      revPE       = 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_revPE_20200804162733_8';  %-maybe wrong alignment
   elseif analysis_id == 34
      image_name  = 'BB_004_3T_cmrr_mbep2d_se_rs_20200804162733_12';
      revPE       = 'BB_004_3T_cmrr_mbep2d_se_revPE_rs_20200804162733_14';
   elseif analysis_id == 35
      image_name  = 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200804162733_15';
      revPE       = 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_revPE_20200804162733_8';  %-maybe wrong alignment
      %%%
   elseif analysis_id == 36
      image_name  = 'BB_005_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200806153736_8';
      revPE       = 'BB_005_3T_cmrr_mbep2d_diff_0pt2_1_revPE_task_20200806153736_14';  %-maybe wrong alignment
   elseif analysis_id == 37
      image_name  = 'BB_005_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200806153736_10';
      revPE       = 'BB_005_3T_cmrr_mbep2d_diff_0pt2_1_revPE_task_20200806153736_14';  %-different b-values between lin and rev! and maybe wrong alignment
   elseif analysis_id == 38
      image_name  = 'BB_005_3T_cmrr_mbep2d_diff_0pt2_1_rs_20200806153736_11';  %-no revPE for resting state
   elseif analysis_id == 39
      image_name  = 'BB_005_3T_cmrr_mbep2d_se_rs_20200806153736_13';  %-no revPE for resting state
      %%%
   elseif analysis_id == 40
      image_name  = 'BB_006_7T_cmrr_mbep2d_diff_b0pt2_1_task_m_seq_20200817103452_19';
      revPE       = 'BB_006_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_20200817103452_8';  %-in fact b-values 0/1!
   elseif analysis_id == 41
      image_name  = 'BB_006_7T_cmrr_mbep2d_se_task_m_seq_20200817103452_21';
      revPE       = 'BB_006_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_20200817103452_8';
   elseif analysis_id == 42
      image_name  = 'BB_006_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20200817103452_23';
      revPE       = 'BB_006_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200817103452_13';  %-different scaling...
   elseif analysis_id == 43
      image_name  = 'BB_006_7T_cmrr_mbep2d_se_rs_20200817103452_25';
      revPE       = 'BB_006_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200817103452_13';
      %%%
   elseif analysis_id == 44
      image_name  = 'BB_006_3T_cmrr_mbep2d_diff_0pt2_1_task_m_seq_with_b_val_0_20200818140942_9';
      revPE       = 'BB_006_3T_cmrr_mbep2d_diff_0pt2_1_revPE_task_b_val_0_20200818140942_15';
      No_dummies_out = 1;
   elseif analysis_id == 45
      image_name  = 'BB_006_3T_cmrr_mbep2d_se_task_m_seq_20200818140942_11';
      revPE       = 'BB_006_3T_cmrr_mbep2d_diff_0pt2_1_revPE_task_b_val_0_20200818140942_15';
   elseif analysis_id == 46
      image_name  = 'BB_006_3T_cmrr_mbep2d_diff_0pt2_1_rs_with_b_val_0_20200818140942_12';
      revPE       = 'BB_006_3T_cmrr_mbep2d_diff_0pt2_1_revPE_rs_b_val_0_20200818140942_16';
      No_dummies_out = 1;
   elseif analysis_id == 47
      image_name  = 'BB_006_3T_cmrr_mbep2d_se_rs_20200818140942_14';
      revPE       = 'BB_006_3T_cmrr_mbep2d_diff_0pt2_1_revPE_rs_b_val_0_20200818140942_16';
      %%%
   elseif analysis_id == 48
      image_name  = 'BB_007_3T_cmrr_mbep2d_se_task_m_seq_20200819141508_9';  %-no proper revPE... other anat image!
   elseif analysis_id == 49
      image_name  = 'BB_007_3T_cmrr_mbep2d_se_rs_20200819141508_14';  %-no proper revPE... other anat image!
      %%%
   elseif analysis_id == 50
      image_name  = 'BB_007_7T_cmrr_mbep2d_diff_b0pt2_1_task_m_seq_20200902101725_24';
      revPE       = 'BB_007_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_20200902101725_8';
   elseif analysis_id == 51
      image_name  = 'BB_007_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20200902101725_26';
      revPE       = 'BB_007_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200902101725_22';  %-different spaces! maybe the subject was taken out and in?
   elseif analysis_id == 52
      image_name  = 'BB_007_7T_cmrr_mbep2d_se_task_m_seq_20200902101725_28';
      revPE       = 'BB_007_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_20200902101725_8';
   elseif analysis_id == 53
      image_name  = 'BB_007_7T_cmrr_mbep2d_se_rs_20200902101725_30';
      revPE       = 'BB_007_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200902101725_22';  %-different spaces! maybe the subject was taken out and in?
      %%%
   elseif analysis_id == 54
      image_name  = 'BB_002_3T_cmrr_mbep2d_se_rs_20200902132205_10';
      revPE       = 'BB_002_3T_cmrr_mbep2d_se_revPE_rs_20200902132205_18';
   elseif analysis_id == 55
      image_name  = 'BB_002_3T_cmrr_mbep2d_diff_0pt2_1_rs_20200902132205_8';
      revPE       = 'BB_002_3T_cmrr_mbep2d_diff_0pt2_1_revPE_rs_20200902132205_19';  %-different scales...
   elseif analysis_id == 56
      image_name  = 'BB_002_3T_cmrr_mbep2d_diff_0pt29_1pt45_rs_20200902132205_15';
      revPE       = 'BB_002_3T_cmrr_mbep2d_diff_0pt29_1pt45_revPE_rs_20200902132205_16';
   elseif analysis_id == 57
      image_name  = 'BB_002_3T_cmrr_mbep2d_diff_0pt2_1_task_m_seq_20200902132205_22';
      revPE       = 'BB_002_3T_cmrr_mbep2d_diff_0pt2_1_revPE_task_20200902132205_23';
      %%%
   elseif analysis_id == 58
      image_name  = 'BB_008_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200910095831_16';
      revPE       = 'BB_008_7T_cmrr_mbep2d_diff_b0_0pt2_1_revPE_task_20200910095831_8';
   elseif analysis_id == 59
      image_name  = 'BB_008_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20200910095831_20';
      revPE       = 'BB_008_7T_cmrr_mbep2d_diff_b0_0pt2_1_revPE_rs_20200910095831_14';
   elseif analysis_id == 60
      image_name  = 'BB_008_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200910095831_18';  %-I adjusted the number of slices removing 4 of them
      revPE       = 'BB_008_7T_cmrr_mbep2d_diff_b0_0pt2_1_revPE_task_20200910095831_8';
   elseif analysis_id == 61
      image_name  = 'BB_008_7T_cmrr_mbep2d_se_rs_20200910095831_22';  %-I adjusted the number of slices removing 4 of them
      revPE       = 'BB_008_7T_cmrr_mbep2d_diff_b0_0pt2_1_revPE_rs_20200910095831_14';
      %%%
   elseif analysis_id == 62
      image_name  = 'BB_008_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200910143135_10';  %-there was poor registration to MNI space when using the T1 map as the brain mask was poor, so I chose a different anat map
      revPE       = 'BB_008_3T_cmrr_mbep2d_diff_0_0pt2_1_task_revPE_20200910143135_11';  %-linPE and revPE are not so similar, as diff=0 is bipolar, SE not
   elseif analysis_id == 63
      image_name  = 'BB_008_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200910143135_8';
      revPE       = 'BB_008_3T_cmrr_mbep2d_diff_0_0pt2_1_task_revPE_20200910143135_11';  %-there was no 0.2... 1x0+6x1000
   elseif analysis_id == 64
      image_name  = 'BB_008_3T_cmrr_mbep2d_diff_0pt2_1_rs_20200910143135_23';
      revPE       = 'BB_008_3T_cmrr_mbep2d_diff_0_0pt2_1_rs_revPE_20200910143135_24';  %-there was no 0.2... 1x0+6x1000
      %%%
   elseif analysis_id == 65
      image_name  = 'BB_009_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200916093315_16';
      revPE       = 'BB_009_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_20200916093315_10';
   elseif analysis_id == 66
      image_name  = 'BB_009_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200916093315_18';
      revPE       = 'BB_009_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_20200916093315_10';  %-different b-values between lin and rev!
   elseif analysis_id == 67
      image_name  = 'BB_009_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20200916093315_20';
      revPE       = 'BB_009_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200916093315_14';
   elseif analysis_id == 68
      image_name  = 'BB_009_7T_cmrr_mbep2d_se_rs_20200916093315_22';
      revPE       = 'BB_009_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_rs_20200916093315_14';  %-different b-values between lin and rev!
      %%%
   elseif analysis_id == 69
      image_name  = 'BB_009_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200918151529_10';
      revPE       = 'BB_009_3T_cmrr_mbep2d_se_revPE_20200918151529_14';
   elseif analysis_id == 70
      image_name  = 'BB_009_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200918151529_8';
      revPE       = 'BB_009_3T_cmrr_mbep2d_diff_0pt2_1_revPE_task_20200918151529_11';
   elseif analysis_id == 71
      image_name  = 'BB_009_3T_cmrr_mbep2d_diff_0pt2_1_rs_20200918151529_15';
      revPE       = 'BB_009_3T_cmrr_mbep2d_diff_0pt2_1_revPE_rs_20200918151529_16';
      %%%
   elseif analysis_id == 72
      image_name  = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200929091524_27';  %-maybe asleep... in zfstat1 one weird slice with a lot of sig
      revPE       = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_MB2_interleaved_bipolar_0pt2';
   elseif analysis_id == 73
      image_name  = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_20200929091524_31';  %-maybe asleep... in zfstat1 some weird clusters
      revPE       = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_MB2_interleaved_bipolar_0pt2';
   elseif analysis_id == 74
      image_name  = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_task_40s_off_40s_on_20200929091524_33';  %-maybe asleep...
      revPE       = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_MB2_interleaved_bipolar_0pt2';
   elseif analysis_id == 75
      image_name  = 'BB_010_7T_cmrr_mbep2d_se_task_18s_off_12s_on_20200929091524_29';  %-asleep...
      revPE       = 'BB_010_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_task_MB2_interleaved_bipolar_0pt2';
      %%%
   elseif analysis_id == 76
      image_name  = 'BB_011_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201001092113_19';
      revPE       = ''; %'BB_011_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201001092113_15';  %-only for b-value 1! very poor SNR! very poor topup output
   elseif analysis_id == 77
      image_name  = 'BB_011_7T_cmrr_mbep2d_diff_b0pt2_1_task_40s_off_40s_on_20201001092113_21';
   elseif analysis_id == 78
      image_name  = 'BB_011_7T_cmrr_mbep2d_diff_b0pt2_1_task_40s_off_40s_on_monopolar_20201001092113_23';
   elseif analysis_id == 79
      image_name  = 'BB_011_7T_cmrr_mbep2d_se_task_18s_off_12s_on_200_vols_20201001092113_25';
   elseif analysis_id == 80
      image_name  = 'BB_011_7T_cmrr_mbep2d_diff_b0pt2_1_rs_300_vols_20201001092113_29';
      revPE       = ''; %'BB_011_7T_cmrr_mbep2d_diff_b0pt2_rs_revPE_20201001092113_17';  %-only for b-value 1! very poor SNR! very poor topup output
   elseif analysis_id == 81
      image_name  = 'BB_011_7T_cmrr_mbep2d_se_rs_300_vols_20201001092113_27';
      %%%
   elseif analysis_id == 82
      image_name  = 'BB_011_3T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201001163631_9';  %-linPE and revPE are very similar, but SDC output looks very good
      revPE       = 'BB_011_3T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201001163631_12';
   elseif analysis_id == 83
      image_name  = 'BB_011_3T_cmrr_mbep2d_diff_b0pt2_1_task_40s_off_40s_on_20201001163631_10';
      revPE       = 'BB_011_3T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201001163631_12';
      %%%
   elseif analysis_id == 84
      image_name  = 'BB_013_3T_cmrr_mbep2d_diff_2pt5iso_b0pt2_1_1206_task_18s_off_12s_on_20201023162033_11_combined';
      revPE       = 'BB_013_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201023162033_9_short';
      No_dummies_out = 6;
   elseif analysis_id == 85
      image_name  = 'BB_013_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20201023162033_14';
      revPE       = 'BB_013_3T_cmrr_mbep2d_diff_rs_2iso_8_revPE_20201023162033_13_short';
   elseif analysis_id == 86
      image_name  = 'BB_013_3T_cmrr_mbep2d_se_rs_20201023162033_15';
      revPE       = 'BB_013_3T_cmrr_mbep2d_diff_rs_2iso_8_revPE_20201023162033_13_short';
      %%%
   elseif analysis_id == 87
      image_name  = 'BB_013_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201022092655_16_combined';
      revPE       = 'BB_013_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201022092655_8_short';
      No_dummies_out = 6;
   elseif analysis_id == 88
      image_name  = 'BB_013_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20201022092655_18';
      revPE       = 'BB_013_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_20201022092655_12_short';
   elseif analysis_id == 89
      image_name  = 'BB_013_7T_cmrr_mbep2d_se_rs_20201022092655_20';
      revPE       = 'BB_013_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_20201022092655_12_short';
      %%%
   elseif analysis_id == 90
      image_name  = 'BB_014_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201029092350_14_combined';
      revPE       = 'BB_014_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201029092350_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 91
      image_name  = 'BB_014_7T_cmrr_mbep2d_se_rs_20201029092350_18';
      revPE       = 'BB_014_7T_cmrr_mbep2d_se_rs_revPE_20201029092350_17';
   elseif analysis_id == 92
      image_name  = 'BB_014_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20201029092350_16';
      revPE       = 'BB_014_7T_cmrr_mbep2d_se_rs_revPE_20201029092350_17';
      %%%
   elseif analysis_id == 93
      image_name  = 'BB_014_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_SENSE1_20201106142023_8_combined';
      revPE       = 'BB_014_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_SENSE1_20201106142023_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 94
      image_name  = 'BB_014_3T_cmrr_mbep2d_se_rs_2iso_SENSE1_20201106142023_12';
      revPE       = 'BB_014_3T_cmrr_mbep2d_se_rs_revPE_20201106142023_10';
   elseif analysis_id == 95
      image_name  = 'BB_014_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20201106142023_11';
      revPE       = 'BB_014_3T_cmrr_mbep2d_se_rs_revPE_20201106142023_10';
      %%%
   elseif analysis_id == 96
      image_name  = 'BB_015_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201103150507_15_combined';
      revPE       = 'BB_015_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201103150507_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 97
      image_name  = 'BB_015_7T_cmrr_mbep2d_se_rs_20201103150507_17';
      revPE       = 'BB_015_7T_cmrr_mbep2d_se_rs_revPE_20201103150507_19';
   elseif analysis_id == 98
      image_name  = 'BB_015_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20201103150507_18';
      revPE       = 'BB_015_7T_cmrr_mbep2d_se_rs_revPE_20201103150507_19';
      %%%
   elseif analysis_id == 99
      image_name  = 'BB_015_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_20201029143232_8_combined';
      revPE       = 'BB_015_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201029143232_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 100
      image_name  = 'BB_015_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20201029143232_11';
      revPE       = 'BB_015_3T_cmrr_mbep2d_diff_rs_2iso_8_revPE_20201029143232_10_short';
   elseif analysis_id == 101
      image_name  = 'BB_015_3T_cmrr_mbep2d_se_rs_20201029143232_12';
      revPE       = 'BB_015_3T_cmrr_mbep2d_diff_rs_2iso_8_revPE_20201029143232_10_short';
      %%%
   elseif analysis_id == 102
      image_name  = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201030143609_9_combined';
      revPE       = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201030143609_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 103
      image_name  = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_20201030143609_13_combined';
      revPE       = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201030143609_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 104
      image_name  = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20201030143609_11';
      revPE       = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_20201030143609_8_short';
   elseif analysis_id == 105
      image_name  = 'BB_016_7T_cmrr_mbep2d_se_rs_20201030143609_12';
      revPE       = 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_20201030143609_8_short';
      %%%
   elseif analysis_id == 106
      image_name  = 'BB_016_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_20201029160019_15_combined';
      revPE       = 'BB_016_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201029160019_17_short';
      No_dummies_out = 6;
   elseif analysis_id == 107
      image_name  = 'BB_016_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20201029160019_21';
      revPE       = 'BB_016_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201029160019_17_short';
   elseif analysis_id == 108
      image_name  = 'BB_016_3T_cmrr_mbep2d_se_rs_20201029160019_20';
      revPE       = 'BB_016_3T_cmrr_mbep2d_se_rs_revPE_20201029160019_19';
   elseif analysis_id == 109
      image_name  = 'BB_016_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20201029160019_18';
      revPE       = 'BB_016_3T_cmrr_mbep2d_se_rs_revPE_20201029160019_19';
      %%%
   elseif analysis_id == 110
      image_name  = 'BB_017_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201106093432_10_combined';
      revPE       = 'BB_017_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201106093432_8_short';
      No_dummies_out = 6;
   elseif analysis_id == 111
      image_name  = 'BB_017_7T_cmrr_mbep2d_se_rs_SENSE1_20201106093432_13';
      revPE       = 'BB_017_7T_cmrr_mbep2d_se_rs_SENSE1_revPE_20201106093432_14';
   elseif analysis_id == 112
      image_name  = 'BB_017_7T_cmrr_mbep2d_diff_b0pt2_1_rs_SENSE1_20201106093432_7';
      revPE       = 'BB_017_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_SENSE1_20201106093432_9_short';
      %%%
   elseif analysis_id == 113
      image_name  = 'BB_001_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_600_18s_off_12s_on_MB2_20201106154548_16';
      revPE       = 'BB_001_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201106154548_8_short';
      No_dummies_out = 1;
   elseif analysis_id == 114
      image_name  = 'BB_001_3T_cmrr_mbep2d_diff_task_2pt5iso_b0_1_600_18s_off_12s_on_20201106154548_15';
      revPE       = 'BB_001_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201106154548_8_short';
      No_dummies_out = 1;
   elseif analysis_id == 115
      image_name  = 'BB_001_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_600_18s_off_12s_on_MB1_8sl_20201106154548_17';
      revPE       = 'BB_001_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201106154548_8_short_cut_to_8sl';
      No_dummies_out = 1;
   elseif analysis_id == 116
      image_name  = 'BB_001_3T_cmrr_mbep2d_se_rs_2iso_20201106154548_10';
      revPE       = 'BB_001_3T_cmrr_mbep2d_se_rs_revPE_20201106154548_11';
   elseif analysis_id == 117
      image_name  = 'BB_001_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20201106154548_12';
      revPE       = 'BB_001_3T_cmrr_mbep2d_se_rs_revPE_20201106154548_11';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 118
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_20201116112841_14';
   elseif analysis_id == 119
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_50_20201116112841_29';
   elseif analysis_id == 120
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_50_20201116112841_30';
   elseif analysis_id == 121
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_60_20201116112841_15';
   elseif analysis_id == 122
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_60_20201116112841_31';
   elseif analysis_id == 123
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_60_int_motion_20201116112841_26';
   elseif analysis_id == 124
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_70_20201116112841_32';
   elseif analysis_id == 125
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_80_20201116112841_33';
   elseif analysis_id == 126
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_90_20201116112841_34';
   elseif analysis_id == 127
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_90_MBHCP_20201116112841_35';
   elseif analysis_id == 128
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_90_MBHCP_TR_2_20201116112841_36';
   elseif analysis_id == 129
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_90_MBHCP_TR_4_20201116112841_37';
   elseif analysis_id == 130
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_FA_90_MBHCP_TR_8_20201116112841_38';
   elseif analysis_id == 131
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_PulseDurationx2_20201116112841_16';
   elseif analysis_id == 132
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_PulseDurationx2_int_m_20201116112841_27';
   elseif analysis_id == 133
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_SENSE1_20201116112841_17';
   elseif analysis_id == 134
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_SENSE1_int_motion_20201116112841_28';
   elseif analysis_id == 135
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_int_motion_20201116112841_25';
   elseif analysis_id == 136
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_130_20201116112841_9';
   elseif analysis_id == 137
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_130_int_motion_20201116112841_20';
   elseif analysis_id == 138
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_60_20201116112841_8';
   elseif analysis_id == 139
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_60_int_motion_20201116112841_19';
   elseif analysis_id == 140
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_20201116112841_7';
   elseif analysis_id == 141
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_PulseDurationx2_20201116112841_13';
   elseif analysis_id == 142
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_PulseDurationx2_int_motion_20201116112841_24';
   elseif analysis_id == 143
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_disabled_grad_rev_20201116112841_10';
   elseif analysis_id == 144
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_disabled_grad_rev_int_motion_20201116112841_21';
   elseif analysis_id == 145
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_disabled_grad_rev_no_fat_supp_int_motion_20201116112841_22';
   elseif analysis_id == 146
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_disabled_grad_rev_no_fat_suppr_20201116112841_11';
   elseif analysis_id == 147
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_int_motion_20201116112841_18';
   elseif analysis_id == 148
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_with_SENSE1_20201116112841_12';
   elseif analysis_id == 149
      image_name  = 'BB_016_7T_PILOT_cmrr_mbep2d_se_rs_FA_90_with_SENSE1_int_motion_20201116112841_23';
      %%%
   elseif analysis_id == 150
      image_name  = 'BB_017_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_20201120153347_10_combined';
      revPE       = 'BB_017_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_20201120153347_9_short';
      No_dummies_out = 8;
   elseif analysis_id == 151
      image_name  = 'BB_017_3T_cmrr_mbep2d_se_rs_2iso_20201120153347_12';
      revPE       = 'BB_017_3T_cmrr_mbep2d_se_rs_revPE_with_SENSE1_20201120153347_7';
   elseif analysis_id == 152
      image_name  = 'BB_017_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20201120153347_13';
      revPE       = 'BB_017_3T_cmrr_mbep2d_se_rs_revPE_with_SENSE1_20201120153347_7';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 153
      image_name  = 'BB_018_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_20210122152445_9_combined';
      revPE       = 'BB_018_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20210122152445_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 154
      image_name  = 'BB_018_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_20210122152445_11_combined';
      revPE       = 'BB_018_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20210122152445_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 155
      image_name  = 'BB_018_7T_cmrr_mbep2d_diff_b0pt2_1_rs_20210122152445_13_short';
      revPE       = 'BB_018_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_bipolar_20210122152445_8_short';
      %%%
   elseif analysis_id == 156
      image_name  = 'BB_018_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_bipolar_20210121150934_7_combined';
      revPE       = 'BB_018_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_bipolar_revPE_20210121150934_14_short';
      No_dummies_out = 8;
   elseif analysis_id == 157
      image_name  = 'BB_018_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_monopolar_20210121150934_9_combined';
      revPE       = 'BB_018_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_bipolar_revPE_20210121150934_14_short';
      No_dummies_out = 8;
   elseif analysis_id == 158
      image_name  = 'BB_018_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_20210121150934_13';  %-monopolar or bipolar? not sure...
      revPE       = 'BB_018_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_revPE_20210121150934_15_short';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 159
      image_name  = 'BB_019_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_20210127143328_12_combined';  %-revPE is AP...
      No_dummies_out = 6;
   elseif analysis_id == 160
      image_name  = 'BB_019_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_2_20210127143328_15_combined';  %-revPE is AP...
      No_dummies_out = 6;
   elseif analysis_id == 161
      image_name  = 'BB_019_7T_cmrr_mbep2d_diff_b0pt2_1_rs_monopolar_20210127143328_11';
      revPE       = 'BB_019_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_monopolar_20210127143328_7_short';
      %%%
   elseif analysis_id == 162
      image_name  = 'BB_019_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_monopolar_20210128131330_7_combined';
      revPE       = 'BB_019_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_monopolar_20210128131330_8_short';
      No_dummies_out = 8;
   elseif analysis_id == 163
      image_name  = 'BB_019_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_bipolar_20210128131330_12_combined';
      revPE       = 'BB_019_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_bipolar_20210128131330_13_short';
      No_dummies_out = 8;
   elseif analysis_id == 164
      image_name  = 'BB_019_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_20210128131330_10';
      revPE       = 'BB_019_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_monopolar_revPE_20210128131330_11_short';  %-wrong scaling used (1000 instead of 200)... only b-values 0 and 1, topup will work on b-values 0 and 0.2...
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 165
      image_name  = 'BB_021_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_monopolar_20210212141923_7_combined';
      revPE       = 'BB_021_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_revPE_monopolar_20210212141923_14_short';
      No_dummies_out = 8;
   elseif analysis_id == 166
      image_name  = 'BB_021_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_bipolar_20210212141923_16_combined';
      revPE       = 'BB_021_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_revPE_monopolar_20210212141923_14_short';
      No_dummies_out = 8;
   elseif analysis_id == 167
      image_name  = 'BB_021_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_20210212141923_12';
      revPE       = 'BB_021_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_revPE_monopolar_20210212141923_17_short';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 168
      image_name  = 'BB_022_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_monopolar_20210212161146_7_combined';
      revPE       = 'BB_022_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_revPE_monopolar_20210212161146_9_short';
      No_dummies_out = 8;
   elseif analysis_id == 169
      image_name  = 'BB_022_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_bipolar_20210212161146_12_combined';
      revPE       = 'BB_022_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_revPE_bipolar_20210212161146_14_short';
      No_dummies_out = 8;
   elseif analysis_id == 170
      image_name  = 'BB_022_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_20210212161146_10';
      revPE       = 'BB_022_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_revPE_20210212161146_11_short';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 171
      image_name  = 'BB_022_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_20210217143129_11_combined';
      revPE       = 'BB_022_7T_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_monopolar_revPE_20210217143129_7_short';
      No_dummies_out = 6;
   elseif analysis_id == 172
      image_name  = 'BB_022_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_20210217143129_14_combined';
      revPE       = 'BB_022_7T_cmrr_mbep2d_diff_b0pt2_task_bipolar_revPE_20210217143129_10_short';
      No_dummies_out = 6;
   elseif analysis_id == 173
      image_name  = 'BB_022_7T_cmrr_mbep2d_diff_b0pt2_1_rs_monopolar_20210217143129_13';
      revPE       = 'BB_022_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_monopolar_20210217143129_9_short';
      %%%
   elseif analysis_id == 174
      image_name  = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_20210219110948_10_combined';
      revPE       = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210219110948_7_short';  %-weird encoding
      No_dummies_out = 6;
   elseif analysis_id == 175
      image_name  = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_1_rs_monopolar_20210219110948_12';
      revPE       = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_monopolar_20210219110948_8_short';
   elseif analysis_id == 176
      image_name  = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_1_rs_bipolar_20210219110948_13';
      revPE       = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_bipolar_20210219110948_9_short';
      %%%
   elseif analysis_id == 177
      image_name  = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_monopolar_1_20210222163533_8';
      revPE       = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_monopolar_20210222163533_9_short';
      No_dummies_out = 1;
   elseif analysis_id == 178
      image_name  = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_604_18s_off_12s_on_monopolar_2_20210222163533_14';
      revPE       = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_monopolar_20210222163533_9_short';
      No_dummies_out = 1;
   elseif analysis_id == 179
      image_name  = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_602_18s_off_12s_on_bipolar_1_20210222163533_11';
      revPE       = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_bipolar_20210222163533_10_short';
      No_dummies_out = 1;
   elseif analysis_id == 180
      image_name  = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_602_18s_off_12s_on_bipolar_2_20210222163533_16';
      revPE       = 'BB_023_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_bipolar_20210222163533_10_short';
      No_dummies_out = 1;
   elseif analysis_id == 181
      image_name  = 'BB_023_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_20210222163533_12';
      revPE       = 'BB_023_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_revPE_20210222163533_15_short';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 182
      image_name  = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_18s_off_12s_on_monopolar_1_20210223162740_7';
      revPE       = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_monopolar_20210223162740_8_short';
      No_dummies_out = 1;
   elseif analysis_id == 183
      image_name  = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_18s_off_12s_on_monopolar_2_20210223162740_13';
      revPE       = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_monopolar_20210223162740_8_short';
      No_dummies_out = 1;
   elseif analysis_id == 184
      image_name  = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_18s_off_12s_on_bipolar_1_20210223162740_9';
      revPE       = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_bipolar_20210223162740_10_short';
      No_dummies_out = 1;
   elseif analysis_id == 185
      image_name  = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_1_18s_off_12s_on_bipolar_2_20210223162740_14';
      revPE       = 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_bipolar_20210223162740_10_short';
      No_dummies_out = 1;
   elseif analysis_id == 186
      image_name  = 'BB_024_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_20210223162740_11';
      revPE       = 'BB_024_3T_cmrr_mbep2d_diff_rs_2iso_b0pt2_1_602_monopolar_revPE_20210223162740_12_short';
      No_dummies_out = 1;
      %%%
   elseif analysis_id == 187
      image_name  = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_1_20210226113202_11';
      revPE       = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210226113202_8_short';
   elseif analysis_id == 188
      image_name  = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_1_20210226113202_12';
      revPE       = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_bipolar_20210226113202_9_short';
   elseif analysis_id == 189
      image_name  = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_1_rs_monopolar_1_20210226113202_10';
      revPE       = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_1_rs_revPE_monopolar_1_20210226113202_7_short';
   elseif analysis_id == 190
      image_name  = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_2_20210226113202_13';
      revPE       = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210226113202_8_short';
   elseif analysis_id == 191
      image_name  = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_2_20210226113202_14';
      revPE       = 'BB_023_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_bipolar_20210226113202_9_short';
      %%%
   elseif analysis_id == 192
      image_name  = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_1_20210226150433_10';
      revPE       = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210226150433_7_short';
   elseif analysis_id == 193
      image_name  = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_1_20210226150433_11';
      revPE       = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_bipolar_20210226150433_8_short';
   elseif analysis_id == 194
      image_name  = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_1_rs_monopolar_20210226150433_12';
      revPE       = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_rs_revPE_monopolar_20210226150433_9_short';
   elseif analysis_id == 195
      image_name  = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_monopolar_2_20210226150433_13';
      revPE       = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210226150433_7_short';
   elseif analysis_id == 196
      image_name  = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_bipolar_2_20210226150433_14';
      revPE       = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_bipolar_20210226150433_8_short';
   end
   
   make_adc       = 0;
   task_analysis  = 0;
   do_topup       = 0;
   
   %-always removing at least two first volumes
   No_dummies_out = No_dummies_out + 2;
   
   if contains(image_name, '_diff_');      make_adc   = 1;     end
   
   if ~isempty(revPE);                     do_topup   = 1;     end
   
   if contains(image_name, 'b0_1');        b_val_diff = 1;     end
   if contains(image_name, '0pt2_1');      b_val_diff = 0.8;   end
   if contains(image_name, '0pt29_1pt45'); b_val_diff = 1.155; end  %-maybe should be 1.16
   
   if contains(image_name, '_task_')
      task_analysis = 1;
      if     contains(image_name, '_18s_off_12s_on_')
         rest_block = 18;
         task_block = 12;
      elseif contains(image_name, '_12s_off_18s_on_')
         rest_block = 12;
         task_block = 18;
      elseif contains(image_name, '_18s_off_2s_on_')
         rest_block = 18;
         task_block = 2;
      elseif contains(image_name, '_40s_off_40s_on_')
         rest_block = 40;
         task_block = 40;
      elseif contains(image_name, '_m_seq_')
         rest_block = 18;
         task_block = 0.1;
      else
         disp('missing paradigm specification!');
         continue
      end
   end
   
   fig_title_text = image_name(1:strfind(image_name, '_2020')-1);
   if isempty(fig_title_text)
      fig_title_text = image_name(1:strfind(image_name, '_2021')-1);
   end
   fig_title_text = strrep(fig_title_text, '_', ' ');
   
   if contains(image_name, '_7T_')
      vol_outliers_thr = vol_outliers_thr*3;
   end
   
   %-identifying the structural (MP2RAGE or MPRAGE) image
   anat_image    = dir(fullfile(path_manage, 'scans', [image_name(1:9) '*mp2rage*']));
   if isempty(anat_image)
      anat_image = dir(fullfile(path_manage, 'scans', [image_name(1:9) '*mprage*']));
   end
   
   anat_name = anat_image(1).name;
   anat_name = anat_name(1:end-7);
   
   if exist(fullfile(path_manage, 'analyses_output', image_name), 'dir') == 7
      disp('directory already exists!!!')
      continue
   end
   
   mkdir(fullfile(path_manage, 'analyses_output', image_name));
   
   copyfile(fullfile(path_manage, 'scans', [image_name '.nii.gz']), fullfile(path_manage, 'analyses_output', image_name, 'dfMRI_raw.nii.gz'));
   
   cd(fullfile(path_manage, 'analyses_output', image_name));
   
   %-open log
   fileID = fopen(fullfile(pwd, 'log.txt'), 'a+');
   fprintf(fileID, [image_name '\n\n']);
   
   %-changing the format to 'FLOAT32', e.g. from 'INT16'
   %-without it, later transformations on the maps are problematic
   system('fslmaths dfMRI_raw dfMRI_raw -odt float');
   
   header_4D  = niftiinfo('dfMRI_raw.nii.gz');
   
   %-due to relaxation status problems in the most outside slices, these slices are being removed
   %-for one image only 4 slices available, so even though the relaxation status problem present, we cannot remove these slices
   %-without these slices (3D) MP-PCA denoising doesnt work
   if ~strcmp(image_name, 'BB_003_3T_cmrr_mbep2d_diff_0pt2_1_TR350_task_18s_off_12s_on_20200710180104_18')
      system(['fslroi dfMRI_raw dfMRI_raw 0 -1 0 -1 1 ' num2str(header_4D.ImageSize(3)-2) ' 0 -1']);
   end
   
   header_4D  = niftiinfo('dfMRI_raw.nii.gz');
   No_slices  = header_4D.ImageSize(3);
   No_volumes = header_4D.ImageSize(4);
   TR         = header_4D.PixelDimensions(4);
   
   system('fslroi dfMRI_raw aux 0 -1 0 -1 0 -1 0 1');
   header_3D  = niftiinfo('aux.nii.gz');
   delete aux.nii.gz
   
   header_3D.Datatype = 'double';
   header_4D.Datatype = 'double';
   
   dfMRI_raw  = niftiread('dfMRI_raw.nii.gz');
   dfMRI_raw  = double(dfMRI_raw);
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-getting rid of the first volumes
   %-1st volume for some scans has to be removed as sometimes it has another b-value than the remaining volumes
   %-(we wanted to cover different b-values - calibration issues)
   
   %-the first three volumes here were for b-value 0, two removed before, one more needed to be removed
   if contains(image_name, 'BB_016_7T_PILOT_cmrr_mbep2d_diff') && No_volumes > 30
      No_dummies_out = No_dummies_out + 2;
   elseif contains(image_name, 'BB_016_7T_PILOT_cmrr_mbep2d_diff')
      No_dummies_out = No_dummies_out + 1;
   end
   
   dfMRI_raw  = dfMRI_raw(:,:,:,No_dummies_out+1:end);
   
   header_4D.ImageSize(4) = header_4D.ImageSize(4) - No_dummies_out;
   niftiwrite(dfMRI_raw, 'dfMRI_raw', header_4D, 'compressed', true);
   
   No_volumes = No_volumes - No_dummies_out;
   
   
   %-identifying volume outliers and removing them
   %-such outliers are primarily a problem at 7T
   %-I am still not sure what causes these outliers (signal drops)
   %-it seems they are correlated with motion
   %-possibly a result of using MB with GRAPPA:
   %-https://practicalfmri.blogspot.com/2012/03/grappa-another-warning-about-motion.html
   %-
   %-at the same time plotting carpet plots showing the distribution of numbers in the images
   
   figure('visible', 'off');
   set(gca,'visible','off');
   set(get(gca,'children'), 'visible', 'off');
   dfMRI_raw_2D          = reshape(dfMRI_raw, size(dfMRI_raw,1)*size(dfMRI_raw,2)*size(dfMRI_raw,3), size(dfMRI_raw, 4));
   mean_dfMRI_raw_2D     = mean(dfMRI_raw_2D);
   save('mean_dfMRI_raw_2D', 'mean_dfMRI_raw_2D');
   
   if make_adc == 0
      
      mean_dfMRI_raw_2D_hpf            = detrend(mean_dfMRI_raw_2D) - mean(detrend(mean_dfMRI_raw_2D)) + mean(mean_dfMRI_raw_2D);
      mean_dfMRI_raw_2D_hpf_w_outliers = mean_dfMRI_raw_2D_hpf;
      
      vol_outliers = find(abs(mean_dfMRI_raw_2D_hpf_w_outliers - median(mean_dfMRI_raw_2D_hpf_w_outliers)) > vol_outliers_thr*median(mean_dfMRI_raw_2D_hpf_w_outliers));
      
      save('vol_outliers', 'vol_outliers');
      
      if ~isempty(vol_outliers)
         mean_dfMRI_raw_2D_hpf(vol_outliers) = NaN;
      end
      
      %-interpolation doesnt work for first and last elements - they cannot be NaN
      non_outliers               = setdiff(1:No_volumes, vol_outliers);
      mean_dfMRI_raw_2D_hpf(1)   = mean_dfMRI_raw_2D_hpf(non_outliers(1));
      mean_dfMRI_raw_2D_hpf(end) = mean_dfMRI_raw_2D_hpf(non_outliers(end));
      dfMRI_raw(:,:,:,1)         = dfMRI_raw(:,:,:,non_outliers(1));
      dfMRI_raw(:,:,:,end)       = dfMRI_raw(:,:,:,non_outliers(end));
      if non_outliers(1) > 1
         non_outliers = [1 non_outliers];
      end
      if non_outliers(end) < No_volumes
         non_outliers = [non_outliers No_volumes];
      end
      for volume_id = 2:No_volumes-1
         if isnan(mean_dfMRI_raw_2D_hpf(volume_id))
            id_before = max(non_outliers(non_outliers < volume_id));
            id_after  = min(non_outliers(non_outliers > volume_id));
            mean_dfMRI_raw_2D_hpf(volume_id) = (mean_dfMRI_raw_2D_hpf(id_before) + mean_dfMRI_raw_2D_hpf(id_after)) / 2;
            dfMRI_raw(:,:,:,volume_id)       = (dfMRI_raw(:,:,:,id_before) + dfMRI_raw(:,:,:,id_after)) / 2;
         end
      end
      
      subplot(2,1,1);
      q_l = quantile(reshape(dfMRI_raw_2D, 1, []), 0.001);
      q_u = quantile(reshape(dfMRI_raw_2D, 1, []), 0.999);
      dfMRI_raw_2D_truncated = dfMRI_raw_2D;
      dfMRI_raw_2D_truncated(dfMRI_raw_2D < q_l) = q_l;
      dfMRI_raw_2D_truncated(dfMRI_raw_2D > q_u) = q_u;
      image(dfMRI_raw_2D_truncated);
      xlim([1 No_volumes]);
      title({fig_title_text; ' '; ' '}, 'FontSize', 8);
      xlabel('volume id');
      ylabel('voxel id');
      
      subplot(2,1,2);
      plot(mean_dfMRI_raw_2D,                'color', [0.8,0.8,0.8]); hold on;
      plot(mean_dfMRI_raw_2D_hpf_w_outliers, 'color', 'b');           hold on;
      plot(mean_dfMRI_raw_2D_hpf,            'color', 'r');
      xlim([1 No_volumes]);
      xlabel('volume id');
      
   else
      
      mean_dfMRI_raw_2D_hpf            = zeros(size(mean_dfMRI_raw_2D));
      mean_dfMRI_raw_2D_hpf(1:2:end)   = detrend(mean_dfMRI_raw_2D(1:2:end)) - mean(detrend(mean_dfMRI_raw_2D(1:2:end))) + mean(mean_dfMRI_raw_2D(1:2:end));
      mean_dfMRI_raw_2D_hpf(2:2:end)   = detrend(mean_dfMRI_raw_2D(2:2:end)) - mean(detrend(mean_dfMRI_raw_2D(2:2:end))) + mean(mean_dfMRI_raw_2D(2:2:end));
      mean_dfMRI_raw_2D_hpf_w_outliers = mean_dfMRI_raw_2D_hpf;
      
      vol_outliers_1 = find(abs(mean_dfMRI_raw_2D_hpf_w_outliers(1:2:end) - median(mean_dfMRI_raw_2D_hpf_w_outliers(1:2:end))) > vol_outliers_thr*median(mean_dfMRI_raw_2D_hpf_w_outliers(1:2:end)));
      vol_outliers_2 = find(abs(mean_dfMRI_raw_2D_hpf_w_outliers(2:2:end) - median(mean_dfMRI_raw_2D_hpf_w_outliers(2:2:end))) > vol_outliers_thr*median(mean_dfMRI_raw_2D_hpf_w_outliers(2:2:end)));
      
      save('vol_outliers_1', 'vol_outliers_1');
      save('vol_outliers_2', 'vol_outliers_2');
      
      if ~isempty(vol_outliers_1)
         mean_dfMRI_raw_2D_hpf(vol_outliers_1*2-1) = NaN;
      end
      if ~isempty(vol_outliers_2)
         mean_dfMRI_raw_2D_hpf(vol_outliers_2*2)   = NaN;
      end
      
      for series = [1 2]
         
         if series == 1
            non_outliers = setdiff(1:No_volumes/2, vol_outliers_1);
            mean_dfMRI_raw_2D_hpf_aux = mean_dfMRI_raw_2D_hpf(1:2:end);
            dfMRI_raw_aux             = dfMRI_raw(:,:,:,1:2:end);
         else
            non_outliers = setdiff(1:No_volumes/2, vol_outliers_2);
            mean_dfMRI_raw_2D_hpf_aux = mean_dfMRI_raw_2D_hpf(2:2:end);
            dfMRI_raw_aux             = dfMRI_raw(:,:,:,2:2:end);
         end
         
         %-interpolation doesnt work for first and last elements - they cannot be NaN
         mean_dfMRI_raw_2D_hpf_aux(1)   = mean_dfMRI_raw_2D_hpf_aux(non_outliers(1));
         mean_dfMRI_raw_2D_hpf_aux(end) = mean_dfMRI_raw_2D_hpf_aux(non_outliers(end));
         dfMRI_raw_aux(:,:,:,1)         = dfMRI_raw_aux(:,:,:,non_outliers(1));
         dfMRI_raw_aux(:,:,:,end)       = dfMRI_raw_aux(:,:,:,non_outliers(end));
         
         if non_outliers(1) > 1
            non_outliers = [1 non_outliers];
         end
         if non_outliers(end) < No_volumes/2
            non_outliers = [non_outliers No_volumes/2];
         end
         
         for volume_id = 2:((No_volumes/2)-1)
            if sum(~isnan(mean_dfMRI_raw_2D_hpf_aux(volume_id))) == 0
               id_before = max(non_outliers(non_outliers < volume_id));
               id_after  = min(non_outliers(non_outliers > volume_id));
               mean_dfMRI_raw_2D_hpf_aux(volume_id) = (mean_dfMRI_raw_2D_hpf_aux(id_before) + mean_dfMRI_raw_2D_hpf_aux(id_after)) / 2;
               dfMRI_raw_aux(:,:,:,volume_id) = (dfMRI_raw_aux(:,:,:,id_before) + dfMRI_raw_aux(:,:,:,id_after)) / 2;
            end
         end
         
         if series == 1
            mean_dfMRI_raw_2D_hpf(1:2:end) = mean_dfMRI_raw_2D_hpf_aux;
            dfMRI_raw(:,:,:,1:2:end)       = dfMRI_raw_aux;
         else
            mean_dfMRI_raw_2D_hpf(2:2:end) = mean_dfMRI_raw_2D_hpf_aux;
            dfMRI_raw(:,:,:,2:2:end)       = dfMRI_raw_aux;
         end
         
      end
      
      subplot(2,2,1);
      q_l = quantile(reshape(dfMRI_raw_2D(:,1:2:end), 1, []), 0.001);
      q_u = quantile(reshape(dfMRI_raw_2D(:,1:2:end), 1, []), 0.999);
      dfMRI_raw_2D_truncated = dfMRI_raw_2D(:,1:2:end);
      dfMRI_raw_2D_truncated(dfMRI_raw_2D_truncated < q_l) = q_l;
      dfMRI_raw_2D_truncated(dfMRI_raw_2D_truncated > q_u) = q_u;
      image(dfMRI_raw_2D_truncated);
      xlim([1 No_volumes/2]);
      title('lower b-value');
      xlabel('volume id');
      ylabel('voxel id');
      set(gca, 'FontSize', 8);
      
      subplot(2,2,2);
      q_l = quantile(reshape(dfMRI_raw_2D(:,2:2:end), 1, []), 0.001);
      q_u = quantile(reshape(dfMRI_raw_2D(:,2:2:end), 1, []), 0.999);
      dfMRI_raw_2D_truncated = dfMRI_raw_2D(:,2:2:end);
      dfMRI_raw_2D_truncated(dfMRI_raw_2D_truncated < q_l) = q_l;
      dfMRI_raw_2D_truncated(dfMRI_raw_2D_truncated > q_u) = q_u;
      image(dfMRI_raw_2D_truncated);
      xlim([1 No_volumes/2]);
      title('higher b-value');
      xlabel('volume id');
      ylabel('voxel id');
      set(gca, 'FontSize', 8);
      
      subplot(2,2,3);
      plot(mean_dfMRI_raw_2D(1:2:end),                'color', [0.8,0.8,0.8]); hold on;
      plot(mean_dfMRI_raw_2D_hpf_w_outliers(1:2:end), 'color', 'b');           hold on;
      plot(mean_dfMRI_raw_2D_hpf(1:2:end),            'color', 'r');
      xlim([1 No_volumes/2]);
      xlabel('volume id');
      set(gca, 'FontSize', 8);
      
      legend('raw signal', 'detrended raw signal', 'detrended raw signal without outliers', 'FontSize', 5, 'Location', 'south');
      legend boxoff
      
      subplot(2,2,4);
      plot(mean_dfMRI_raw_2D(2:2:end),                'color', [0.8,0.8,0.8]); hold on;
      plot(mean_dfMRI_raw_2D_hpf_w_outliers(2:2:end), 'color', 'b');           hold on;
      plot(mean_dfMRI_raw_2D_hpf(2:2:end),            'color', 'r');
      xlim([1 No_volumes/2]);
      xlabel('volume id');
      set(gca, 'FontSize', 8);
      
      sgtitle(fig_title_text, 'FontSize', 8);
      
   end
   
   if analysis_id < 10
      print_to_svg_to_pdf(['dfMRI_raw_carpet_plot_0' num2str(analysis_id)]);
   else
      print_to_svg_to_pdf(['dfMRI_raw_carpet_plot_'  num2str(analysis_id)]);
   end
   
   niftiwrite(dfMRI_raw, 'dfMRI_raw', header_4D, 'compressed', true);
   
   clear dfMRI_raw_2D dfMRI_raw_2D_truncated dfMRI_raw_aux
   
   
   if (make_adc == 0 &&                         length(vol_outliers)   /No_volumes > frac_outliers_discard) || ...
         (make_adc == 1 && (length(vol_outliers_1)+length(vol_outliers_2))/No_volumes > frac_outliers_discard)
      fprintf(fileID, 'too many volume outliers, stopping the analysis of the image');
      continue
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-PCA denoising
   
   %-at the end copying PCA denoising output to a folder dedicated to these files for all subjects
   %-important as I often repeat analyses, and PCA denoising takes much time
   PCA_files_to_copy = {'PCA_denoising_sigma.nii.gz'; 'dfMRI_denoised.nii.gz'; 'PCA_denoising_sigma_part_2.nii.gz'; 'PCA_denoising_npars.nii.gz'; 'PCA_denoising_npars_part_2.nii.gz'};
   
   %-performing PCA denoising only if no previous output available
   if exist(fullfile(path_manage, 'PCA_denoising', image_name, 'dfMRI_denoised.nii.gz'), 'file') ~= 2
      
      for PCA_denoising_part = 1:2
         
         if     make_adc == 0 && PCA_denoising_part == 1
            dfMRI_raw_aux = dfMRI_raw;
         elseif make_adc == 0 && PCA_denoising_part == 2
            continue
         elseif make_adc == 1 && PCA_denoising_part == 1
            dfMRI_raw_aux = dfMRI_raw(:,:,:,1:2:end);
         elseif make_adc == 1 && PCA_denoising_part == 2
            dfMRI_raw_aux = dfMRI_raw(:,:,:,2:2:end);
         end
         
         size_vector          = size(dfMRI_raw_aux);
         size_vector_aux      = size_vector;
         size_vector_aux(1:3) = size_vector_aux(1:3) + 4;
         
         %-denoising in the entire volume, no real brain mask is used
         mask           = ones(size_vector_aux(1:3));
         
         image_aux      = zeros(size_vector_aux);
         
         image_aux(3:end-2,3:end-2,1:2,  :)     = dfMRI_raw_aux(:,:,end-1:end,:);
         image_aux(3:end-2,3:end-2,3:end-2, :)  = dfMRI_raw_aux;
         image_aux(3:end-2,3:end-2,end-1:end,:) = dfMRI_raw_aux(:,:,1:2,:);
         
         image_aux(1:2,:,:,:)       = image_aux(end-3:end-2,:,:,:);
         image_aux(end-1:end,:,:,:) = image_aux(3:4,:,:,:);
         image_aux(:,1:2,:,:)       = image_aux(:,end-3:end-2,:,:);
         image_aux(:,end-1:end,:,:) = image_aux(:,3:4,:,:);
         
         dfMRI_denoised      = zeros(size_vector);
         PCA_denoising_sigma = zeros(size_vector(1:3));
         PCA_denoising_npars = zeros(size_vector(1:3));
         
         %-to use less RAM, I run the denoising on each slice separately
         
         for slice_id = 1:No_slices
            
            [dfMRI_denoised_block, PCA_denoising_sigma_block, PCA_denoising_npars_block] = MPdenoising(image_aux(:,:,slice_id:(slice_id+4),:), mask(:,:,slice_id:(slice_id+4)), [5 5 5], 'full');
            dfMRI_denoised     (:,:,slice_id,:) = dfMRI_denoised_block     (3:(end-2),3:(end-2),3,:);
            PCA_denoising_sigma(:,:,slice_id)   = PCA_denoising_sigma_block(3:(end-2),3:(end-2),3);
            PCA_denoising_npars(:,:,slice_id)   = PCA_denoising_npars_block(3:(end-2),3:(end-2),3);
            
         end
         
         clear dfMRI_denoised_block
         
         if PCA_denoising_part == 1
            dfMRI_denoised_part_1 = dfMRI_denoised;
            niftiwrite(PCA_denoising_sigma, 'PCA_denoising_sigma',        header_3D, 'compressed', true);
            niftiwrite(PCA_denoising_npars, 'PCA_denoising_npars',        header_3D, 'compressed', true);
         elseif PCA_denoising_part == 2
            dfMRI_denoised_combined = zeros(size(dfMRI_raw));
            dfMRI_denoised_combined(:,:,:,1:2:end) = dfMRI_denoised_part_1;
            dfMRI_denoised_combined(:,:,:,2:2:end) = dfMRI_denoised;
            dfMRI_denoised = dfMRI_denoised_combined;
            niftiwrite(PCA_denoising_sigma, 'PCA_denoising_sigma_part_2', header_3D, 'compressed', true);
            niftiwrite(PCA_denoising_npars, 'PCA_denoising_npars_part_2', header_3D, 'compressed', true);
         end
         
      end
      
      clear dfMRI_denoised_part_1 dfMRI_denoised_combined dfMRI_raw_aux
      
      niftiwrite(dfMRI_denoised, 'dfMRI_denoised', header_4D, 'compressed', true);
      
      %-copyfile only works if the destination directory already exists
      mkdir(fullfile(path_manage, 'PCA_denoising', image_name));
      
      for file_id = 1:length(PCA_files_to_copy)
         file_to_copy = PCA_files_to_copy(file_id);
         file_to_copy = file_to_copy{1};
         if exist(file_to_copy, 'file') == 2
            copyfile(file_to_copy, fullfile(path_manage, 'PCA_denoising', image_name, file_to_copy));
         end
      end
      
   else
      
      for file_id = 1:length(PCA_files_to_copy)
         file_to_copy = PCA_files_to_copy(file_id);
         file_to_copy = file_to_copy{1};
         if exist(fullfile(path_manage, 'PCA_denoising', image_name, file_to_copy), 'file') == 2
            copyfile(fullfile(path_manage, 'PCA_denoising', image_name, file_to_copy), file_to_copy);
         end
      end
      
   end
   
   dfMRI_denoised      = niftiread(fullfile(path_manage, 'PCA_denoising', image_name, 'dfMRI_denoised'));
   PCA_denoising_sigma = niftiread(fullfile(path_manage, 'PCA_denoising', image_name, 'PCA_denoising_sigma'));
   PCA_denoising_npars = niftiread(fullfile(path_manage, 'PCA_denoising', image_name, 'PCA_denoising_npars'));
   
   %-saving PCA denoising residuals
   PCA_denoising_residuals = dfMRI_denoised - dfMRI_raw;
   niftiwrite(double(PCA_denoising_residuals), 'PCA_denoising_residuals', header_4D, 'compressed', true);
   
   %-even though residual maps seem like pure noise (when looking volume after volume), the means show spatial inhomogeneities
   system('fslmaths PCA_denoising_residuals -Tmean PCA_denoising_residuals_mean');
   
   clear PCA_denoising_residuals image_aux
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-Gibbs unringing
   dfMRI_denoised_unringed = unring(dfMRI_denoised);
   niftiwrite(dfMRI_denoised_unringed, 'dfMRI_denoised_unringed', header_4D, 'compressed', true);
   clear dfMRI_denoised_unringed
   
   system('fslmaths dfMRI_denoised_unringed      -sub dfMRI_denoised dfMRI_denoised_unringed_diff');
   system('fslmaths dfMRI_denoised_unringed_diff -Tmean              dfMRI_denoised_unringed_diff_mean');
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-calculate tSNR and SNR maps
   
   system('rm -rf SNR');
   mkdir('SNR');
   
   %-for tSNR maps, we need to remove low-frequency trends
   
   dfMRI_raw_hpf         = reshape(dfMRI_raw,      numel(dfMRI_raw)     /size(dfMRI_raw,      4), size(dfMRI_raw,      4))';
   dfMRI_denoised_hpf    = reshape(dfMRI_denoised, numel(dfMRI_denoised)/size(dfMRI_denoised, 4), size(dfMRI_denoised, 4))';
   if make_adc == 0
      dfMRI_raw_hpf      = detrend(dfMRI_raw_hpf)      - mean(detrend(dfMRI_raw_hpf))      + mean(dfMRI_raw_hpf);
      dfMRI_denoised_hpf = detrend(dfMRI_denoised_hpf) - mean(detrend(dfMRI_denoised_hpf)) + mean(dfMRI_denoised_hpf);
   else
      dfMRI_raw_hpf     (1:2:end,:) = detrend(dfMRI_raw_hpf     (1:2:end,:)) - mean(detrend(dfMRI_raw_hpf     (1:2:end,:))) + mean(dfMRI_raw_hpf     (1:2:end,:));
      dfMRI_raw_hpf     (2:2:end,:) = detrend(dfMRI_raw_hpf     (2:2:end,:)) - mean(detrend(dfMRI_raw_hpf     (2:2:end,:))) + mean(dfMRI_raw_hpf     (2:2:end,:));
      dfMRI_denoised_hpf(1:2:end,:) = detrend(dfMRI_denoised_hpf(1:2:end,:)) - mean(detrend(dfMRI_denoised_hpf(1:2:end,:))) + mean(dfMRI_denoised_hpf(1:2:end,:));
      dfMRI_denoised_hpf(2:2:end,:) = detrend(dfMRI_denoised_hpf(2:2:end,:)) - mean(detrend(dfMRI_denoised_hpf(2:2:end,:))) + mean(dfMRI_denoised_hpf(2:2:end,:));
   end
   dfMRI_raw_hpf         = dfMRI_raw_hpf';
   dfMRI_denoised_hpf    = dfMRI_denoised_hpf';
   dfMRI_raw_hpf         = reshape(dfMRI_raw_hpf,      size(dfMRI_raw));
   dfMRI_denoised_hpf    = reshape(dfMRI_denoised_hpf, size(dfMRI_denoised));
   
   niftiwrite(dfMRI_raw_hpf,      'SNR/dfMRI_raw_hpf',      header_4D, 'compressed', true);
   niftiwrite(dfMRI_denoised_hpf, 'SNR/dfMRI_denoised_hpf', header_4D, 'compressed', true);
   
   if make_adc == 0
      
      system('fslmaths dfMRI_raw      -Tmean SNR/dfMRI_raw_mean');
      system('fslmaths dfMRI_denoised -Tmean SNR/dfMRI_denoised_mean');
      
      %-calculate tSNR before PCA denoising
      system( 'fslmaths SNR/dfMRI_raw_hpf       -Tmean                           SNR/dfMRI_raw_mean');
      system( 'fslmaths SNR/dfMRI_raw_hpf       -Tstd                            SNR/dfMRI_raw_std');
      system(['fslmaths SNR/dfMRI_raw_std       -mul ' num2str(sqrt(TR))   '     SNR/dfMRI_raw_std_TR']);
      system( 'fslmaths SNR/dfMRI_raw_mean      -div   SNR/dfMRI_raw_std_TR      SNR/tSNR_before_PCA_denoising');
      
      %-calculate tSNR after PCA denoising
      system( 'fslmaths SNR/dfMRI_denoised_hpf  -Tmean                           SNR/dfMRI_denoised_mean');
      system( 'fslmaths SNR/dfMRI_denoised_hpf  -Tstd                            SNR/dfMRI_denoised_std');
      system(['fslmaths SNR/dfMRI_denoised_std  -mul ' num2str(sqrt(TR))   '     SNR/dfMRI_denoised_std_TR']);
      system( 'fslmaths SNR/dfMRI_denoised_mean -div   SNR/dfMRI_denoised_std_TR SNR/tSNR_after_PCA_denoising');
      
   else
      
      header_4D_aux = header_4D;
      header_4D_aux.ImageSize(4) = floor(header_4D_aux.ImageSize(4)/2);
      
      niftiwrite(double(dfMRI_raw_hpf     (:,:,:,1:2:(No_volumes-1))), 'SNR/dfMRI_raw_hpf_b_0',      header_4D_aux, 'compressed', true);
      niftiwrite(double(dfMRI_raw_hpf     (:,:,:,2:2:No_volumes)),     'SNR/dfMRI_raw_hpf_b_1',      header_4D_aux, 'compressed', true);
      niftiwrite(double(dfMRI_denoised_hpf(:,:,:,1:2:(No_volumes-1))), 'SNR/dfMRI_denoised_hpf_b_0', header_4D_aux, 'compressed', true);
      niftiwrite(double(dfMRI_denoised_hpf(:,:,:,2:2:No_volumes)),     'SNR/dfMRI_denoised_hpf_b_1', header_4D_aux, 'compressed', true);
      
      %-b value 0: calculate tSNR before PCA denoising
      system( 'fslmaths SNR/dfMRI_raw_hpf_b_0       -Tmean                               SNR/dfMRI_raw_mean_b_0');
      system( 'fslmaths SNR/dfMRI_raw_hpf_b_0       -Tstd                                SNR/dfMRI_raw_std_b_0');
      system(['fslmaths SNR/dfMRI_raw_std_b_0       -mul ' num2str(sqrt(TR))      '      SNR/dfMRI_raw_std_TR_b_0']);
      system( 'fslmaths SNR/dfMRI_raw_mean_b_0      -div   SNR/dfMRI_raw_std_TR_b_0      SNR/tSNR_before_PCA_denoising_b_0');
      
      %-b value 0: calculate tSNR after PCA denoising
      system( 'fslmaths SNR/dfMRI_denoised_hpf_b_0  -Tmean                               SNR/dfMRI_denoised_mean_b_0');
      system( 'fslmaths SNR/dfMRI_denoised_hpf_b_0  -Tstd                                SNR/dfMRI_denoised_std_b_0');
      system(['fslmaths SNR/dfMRI_denoised_std_b_0  -mul ' num2str(sqrt(TR))      '      SNR/dfMRI_denoised_std_TR_b_0']);
      system( 'fslmaths SNR/dfMRI_denoised_mean_b_0 -div   SNR/dfMRI_denoised_std_TR_b_0 SNR/tSNR_after_PCA_denoising_b_0');
      
      %-b value 1: calculate tSNR before PCA denoising
      system( 'fslmaths SNR/dfMRI_raw_hpf_b_1       -Tmean                               SNR/dfMRI_raw_mean_b_1');
      system( 'fslmaths SNR/dfMRI_raw_hpf_b_1       -Tstd                                SNR/dfMRI_raw_std_b_1');
      system(['fslmaths SNR/dfMRI_raw_std_b_1       -mul ' num2str(sqrt(TR))      '      SNR/dfMRI_raw_std_TR_b_1']);
      system( 'fslmaths SNR/dfMRI_raw_mean_b_1      -div   SNR/dfMRI_raw_std_TR_b_1      SNR/tSNR_before_PCA_denoising_b_1');
      
      %-b value 1: calculate tSNR after PCA denoising
      system( 'fslmaths SNR/dfMRI_denoised_hpf_b_1  -Tmean                               SNR/dfMRI_denoised_mean_b_1');
      system( 'fslmaths SNR/dfMRI_denoised_hpf_b_1  -Tstd                                SNR/dfMRI_denoised_std_b_1');
      system(['fslmaths SNR/dfMRI_denoised_std_b_1  -mul ' num2str(sqrt(TR))      '      SNR/dfMRI_denoised_std_TR_b_1']);
      system( 'fslmaths SNR/dfMRI_denoised_mean_b_1 -div   SNR/dfMRI_denoised_std_TR_b_1 SNR/tSNR_after_PCA_denoising_b_1');
      
   end
   
   %-calculate SNR before PCA denoising
   %-this is for all volumes, but could be for one volume only
   system('fslmaths dfMRI_raw                -div PCA_denoising_sigma SNR/SNR_before_PCA_denoising');
   
   %-calculate SNR after PCA denoising
   %-sqrt(M/P) increase following Eq. 9 from Veraart ea. 2016
   %-Eq. 9 assumes there are more voxels in the local neighbourhood than there are volumes!
   %-this pipeline works with 5x5x5 local neighbourhoods (125 voxels)
   %-now for diffusion data, PCA denoising runs separately across the two b-values, so 'PCA_SNR_improvement' calculation is wrong!
   if No_volumes < 125
      system(['fslmaths PCA_denoising_npars  -div ' num2str(No_volumes) '     SNR/PCA_SNR_improvement']);
   else
      system( 'fslmaths PCA_denoising_npars  -div 125                         SNR/PCA_SNR_improvement');
   end
   system('fslmaths SNR/PCA_SNR_improvement  -recip                           SNR/PCA_SNR_improvement');
   system('fslmaths SNR/PCA_SNR_improvement  -sqrt                            SNR/PCA_SNR_improvement');
   system('fslmaths SNR/SNR_before_PCA_denoising -mul SNR/PCA_SNR_improvement SNR/SNR_after_PCA_denoising');
   
   close all;
   figure('rend', 'painters', 'pos', [0 0 400 200], 'Visible', 'off');
   subplot(2,1,1);
   plot(reshape(dfMRI_raw(55,55,4,:), 1, [])); hold on;
   plot(reshape(dfMRI_raw_hpf(55,55,4,:), 1, []));
   title('dfMRI raw');
   subplot(2,1,2);
   plot(reshape(dfMRI_denoised(55,55,4,:), 1, []));     hold on;
   plot(reshape(dfMRI_denoised_hpf(55,55,4,:), 1, []));
   title('dfMRI denoised');
   print('hpf', '-dpng');
   
   
   clear dfMRI_raw dfMRI_raw_hpf dfMRI_denoised dfMRI_denoised_hpf
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-susceptibility distortion correction
   
   if do_topup == 1
      
      if task_analysis == 1
         PE_dim = 'LR';
      else
         PE_dim = 'AP';
      end
      
      %-topup exceptions
      if strcmp(   revPE, 'BB_001_7T_cmrr_mbep2d_se_revPE_20200611094346_17') || ...
            strcmp(revPE, 'BB_001_7T_cmrr_mbep2d_diff_b0_1_bipolar_revPE_20200611094346_13') || ...
            strcmp(revPE, 'BB_002_7T_cmrr_mbep2d_diff_b0_1_revPE_20200630150831_14') || ...
            strcmp(revPE, 'BB_002_7T_cmrr_mbep2d_se_revPE_20200630150831_16') || ...
            strcmp(revPE, 'BB_003_7T_cmrr_mbep2d_diff_b0pt2_true_1_revPE_20200702093313_12') || ...
            strcmp(revPE, 'BB_003_7T_cmrr_mbep2d_se_revPE_20200702093313_10') || ...
            strcmp(revPE, 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_20200706151255_8') || ...
            strcmp(revPE, 'BB_004_7T_cmrr_mbep2d_se_revPE_20200706151255_10') || ...
            strcmp(revPE, 'BB_004_7T_cmrr_mbep2d_diff_b0pt29_1pt45_revPE_20200706151255_12') || ...
            strcmp(revPE, 'BB_005_7T_cmrr_mbep2d_diff_b0pt2_1_revPE_20200731112951_3') || ...
            strcmp(revPE, 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_revPE_20200804162733_8')
         PE_dim = 'LR';
      elseif analysis_id < 36
         PE_dim = 'AP';
      elseif strcmp(revPE, 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210219110948_7_short') || ...
            strcmp(revPE, 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_monopolar_20210223162740_8_short') || ...
            strcmp(revPE, 'BB_024_3T_cmrr_mbep2d_diff_task_2pt5iso_b0pt2_8_revPE_bipolar_20210223162740_10_short')
         PE_dim = 'AP';
      end
      
      %-checking if topup already used for this revPE image
      if isempty(dir(fullfile(path_manage, 'topup_output', revPE)))
         
         copyfile(fullfile(path_manage, 'scans', [revPE '.nii.gz']), fullfile(pwd, [revPE '.nii.gz']));
         
         linPE_header    = niftiinfo('dfMRI_raw.nii.gz');
         revPE_header    = niftiinfo([revPE '.nii.gz']);
         
         linPE_No_slices = linPE_header.ImageSize(3);
         revPE_No_slices = revPE_header.ImageSize(3);
         
         linPE_img       = niftiread('dfMRI_raw.nii.gz');
         
         %-it is tricky to use niftiwrite for a 3D map using a header from a 4D map,
         %-so for three revPE images with only one volume, I copied/doubled this volume, so there are two volumes:
         %-BB_008_3T_cmrr_mbep2d_diff_0_0pt2_1_rs_revPE_20200910143135_24.nii.gz
         %-BB_008_3T_cmrr_mbep2d_diff_0_0pt2_1_task_revPE_20200910143135_11.nii.gz
         %-pilot_RS_cmrr_mbep2d_b0_revPE_20200217155721_8.nii.gz
         
         linPE_header.ImageSize(4) = revPE_header.ImageSize(4);
         
         %-for these images, accidentally revPE acquired using the higher b-value only...
         if strcmp(revPE, 'BB_011_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_20201001092113_15') || ...
               strcmp(revPE, 'BB_011_7T_cmrr_mbep2d_diff_b0pt2_rs_revPE_20201001092113_17')
            niftiwrite(linPE_img(:,:,:,2:2:(2*revPE_header.ImageSize(4))), 'linPE_vols', linPE_header, 'compressed', true);
         else
            niftiwrite(linPE_img(:,:,:,1:2:(2*revPE_header.ImageSize(4))), 'linPE_vols', linPE_header, 'compressed', true);
         end
         
         clear linPE_img
         
         %-would be better to include motion correction here
         system( 'fslmaths linPE_vols -Tmean linPE_one_vol');
         system(['fslmaths ' revPE  ' -Tmean revPE_one_vol']);
         
         if     linPE_No_slices < revPE_No_slices
            system(['fslroi revPE_one_vol revPE_one_vol 0 -1 0 -1 ' num2str((revPE_No_slices-linPE_No_slices)/2) ' ' num2str(linPE_No_slices) ' 0 -1']);
         elseif linPE_No_slices > revPE_No_slices
            system(['fslroi linPE_one_vol linPE_one_vol 0 -1 0 -1 ' num2str((linPE_No_slices-revPE_No_slices)/2) ' ' num2str(revPE_No_slices) ' 0 -1']);
         end
         
         system( 'fslmerge -t linPE_and_revPE linPE_one_vol revPE_one_vol');
         
         system( 'fslroi linPE_and_revPE linPE_and_revPE_first 0 -1 0 -1 0 1 0 -1');
         
         system(['fslroi linPE_and_revPE linPE_and_revPE_last  0 -1 0 -1 ' num2str(No_slices-1) ' 1 0 -1']);
         
         system( 'fslmerge -z linPE_and_revPE_padded linPE_and_revPE_first linPE_and_revPE linPE_and_revPE_last');
         
         %-https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/blob/master/bb_data/topup_b02b0.cnf
         system(['topup --imain=linPE_and_revPE_padded --datain=' fullfile(path_manage_software, ['acqp_' PE_dim '.txt']) ' --out=topup_results --iout=linPE_and_revPE_tu --config=' fullfile(path_manage_software, 'topup.cnf') ' -v']);
         
         mkdir(fullfile(path_manage, 'topup_output', revPE));
         movefile('topup_results_fieldcoef.nii.gz', fullfile(path_manage, 'topup_output', revPE, 'topup_results_fieldcoef.nii.gz'));
         movefile('topup_results_movpar.txt',       fullfile(path_manage, 'topup_output', revPE, 'topup_results_movpar.txt'));
         
         files = dir('linPE*');
         for file_id = 1:length(files)
            movefile(files(file_id).name, fullfile(path_manage, 'topup_output', revPE, files(file_id).name));
         end
         
      end
      
      topup_results = dir(fullfile(path_manage, 'topup_output', revPE));
      
      system( 'fslroi      dfMRI_denoised_unringed        dfMRI_denoised_unringed_first 0 -1 0 -1 0                        1 0 -1');
      system(['fslroi      dfMRI_denoised_unringed        dfMRI_denoised_unringed_last  0 -1 0 -1 ' num2str(No_slices-1) ' 1 0 -1']);
      system( 'fslmerge -z dfMRI_denoised_unringed_padded dfMRI_denoised_unringed_first dfMRI_denoised_unringed dfMRI_denoised_unringed_last');
      
      system(['applytopup -i dfMRI_denoised_unringed_padded --inindex=1 --datain=' fullfile(path_manage_software, ['acqp_' PE_dim '.txt']) ' --topup=' fullfile(topup_results(1).folder, 'topup_results') ' --out=dfMRI_denoised_unringed_padded_uw --method=jac']);
      
      system(['fslroi dfMRI_denoised_unringed_padded_uw dfMRI_denoised_unringed_sdc 0 -1 0 -1 1 ' num2str(No_slices) ' 0 -1']);
      
   else
      
      %-if no reversed PE image available, susceptibility distortion correction is not performed
      %-however, in the next step an image with the suffix '_sdc' is needed...
      copyfile('dfMRI_denoised_unringed.nii.gz', 'dfMRI_denoised_unringed_sdc.nii.gz')
      
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-realign (motion correction)
   
   system('gunzip dfMRI_denoised_unringed_sdc.nii.gz');
   
   if make_adc == 0
      
      clear jobs;
      jobs{1}.spatial{1}.realign{1}.estwrite.data{1} = cellstr('dfMRI_denoised_unringed_sdc.nii');
      spm_jobman('run', jobs);
      
      system('gzip rdfMRI_denoised_unringed_sdc.nii');
      
   else
      
      header_4D_aux = header_4D;
      header_4D_aux.ImageSize(4) = floor(header_4D_aux.ImageSize(4)/2);
      
      dfMRI_denoised_unringed_sdc     = niftiread('dfMRI_denoised_unringed_sdc.nii');
      dfMRI_denoised_unringed_sdc_b_0 = dfMRI_denoised_unringed_sdc(:,:,:,1:2:(end-1));
      dfMRI_denoised_unringed_sdc_b_1 = dfMRI_denoised_unringed_sdc(:,:,:,2:2:end);
      
      niftiwrite(dfMRI_denoised_unringed_sdc_b_0, 'dfMRI_denoised_unringed_sdc_b_0', header_4D_aux);
      niftiwrite(dfMRI_denoised_unringed_sdc_b_1, 'dfMRI_denoised_unringed_sdc_b_1', header_4D_aux);
      
      clear dfMRI_denoised_unringed_sdc dfMRI_denoised_unringed_sdc_b_0 dfMRI_denoised_unringed_sdc_b_1
      
      clear jobs;
      jobs{1}.spatial{1}.realign{1}.estwrite.data{1} = cellstr('dfMRI_denoised_unringed_sdc_b_0.nii');
      spm_jobman('run', jobs);
      
      clear jobs;
      jobs{1}.spatial{1}.realign{1}.estwrite.data{1} = cellstr('dfMRI_denoised_unringed_sdc_b_1.nii');
      spm_jobman('run', jobs);
      
      system('fslmaths rdfMRI_denoised_unringed_sdc_b_0 -Tmean rdfMRI_denoised_unringed_sdc_b_0_mean');
      system('fslmaths rdfMRI_denoised_unringed_sdc_b_1 -Tmean rdfMRI_denoised_unringed_sdc_b_1_mean');
      
      %-moving the images for the higher b-value to the space of the images for the lower b-value
      %-the output is a bit stochastic - different runs lead to slightly different registrations...
      %-this step should be improved!!!!!!!!!!!!
      
      %-using ANTs for rigid registration
      system(['antsRegistration --dimensionality 3 --float 0 ' ...
         ' --output [rdfMRI_denoised_unringed_sdc_b_1_mean_in_0, rdfMRI_denoised_unringed_sdc_b_1_mean_in_0.nii] ' ...
         ' --interpolation Linear ' ...
         ' --winsorize-image-intensities [0.01,0.99] ' ...
         ' --use-histogram-matching 0 ' ...
         ' --transform Rigid[0.0001] ' ...
         ' --metric MI[rdfMRI_denoised_unringed_sdc_b_0_mean.nii.gz,rdfMRI_denoised_unringed_sdc_b_1_mean.nii.gz,1,32,Regular,0.25] ' ...
         ' --convergence [10000x10000x10000x10000,1e-6,10] ' ...
         ' --shrink-factors 4x4x2x1 ' ...
         ' --smoothing-sigmas 2x2x1x0vox']);
      
      system(['antsApplyTransforms -r rdfMRI_denoised_unringed_sdc_b_0_mean.nii.gz' ...
         ' -i rdfMRI_denoised_unringed_sdc_b_1.nii -e 3 -t rdfMRI_denoised_unringed_sdc_b_1_mean_in_00GenericAffine.mat -v -o rdfMRI_denoised_unringed_sdc_b_1.nii']);
      
      rdfMRI_denoised_unringed_sdc = zeros(header_4D.ImageSize);
      rdfMRI_denoised_unringed_sdc(:,:,:,1:2:(end-1)) = niftiread('rdfMRI_denoised_unringed_sdc_b_0');
      rdfMRI_denoised_unringed_sdc(:,:,:,2:2:end)     = niftiread('rdfMRI_denoised_unringed_sdc_b_1');
      
      niftiwrite(rdfMRI_denoised_unringed_sdc, 'rdfMRI_denoised_unringed_sdc', header_4D, 'compressed', true);
      
   end
   
   %-replacing NaNs (resulting from 'realign') with zeros
   system('fslmaths rdfMRI_denoised_unringed_sdc -nan rdfMRI_denoised_unringed_sdc');
   system('gunzip   rdfMRI_denoised_unringed_sdc.nii.gz');
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-ideally, slice timing correction should be performed here, but:
   %-(1) for our data there are MB/GRAPPA artifacts -> some slices/volumes have much lower signal
   %-(2) because of SAR problems, different numbers of slices, not always 16
   
   if make_adc == 0
      
      copyfile('rdfMRI_denoised_unringed_sdc.nii', 'ardfMRI_denoised_unringed_sdc.nii');
      
      system('gzip ardfMRI_denoised_unringed_sdc.nii');
      
   else
      
      %-for later ADC analyses, merging the images pair-wise, so that two b-values in one volume
      
      rdfMRI_denoised_unringed_sdc = niftiread('rdfMRI_denoised_unringed_sdc.nii');
      size_aux                     = size(rdfMRI_denoised_unringed_sdc);
      size_aux(3)                  = 2*size_aux(3);
      %-sometimes possible the number of volumes is odd
      size_aux(4)                  = floor(size_aux(4)/2);
      rdfMRI_denoised_unringed_sdc_paired                        = zeros(size_aux);
      rdfMRI_denoised_unringed_sdc_paired(:,:,1:No_slices,    :) = rdfMRI_denoised_unringed_sdc(:,:,:,1:2:(No_volumes-1));
      rdfMRI_denoised_unringed_sdc_paired(:,:,No_slices+1:end,:) = rdfMRI_denoised_unringed_sdc(:,:,:,2:2:No_volumes);
      header_4D_aux                = header_4D;
      header_4D_aux.ImageSize(3:4) = size_aux(3:4);
      niftiwrite(rdfMRI_denoised_unringed_sdc_paired, 'rdfMRI_denoised_unringed_sdc_paired', header_4D_aux);
      
      copyfile('rdfMRI_denoised_unringed_sdc_paired.nii', 'ardfMRI_denoised_unringed_sdc_paired.nii');
      
      ardfMRI_denoised_unringed_sdc_paired       = niftiread('ardfMRI_denoised_unringed_sdc_paired');
      header_4D.ImageSize(4)                     = floor(No_volumes/2);
      ardfMRI_denoised_unringed_sdc_b_0          = zeros(header_4D.ImageSize);
      ardfMRI_denoised_unringed_sdc_b_1          = zeros(header_4D.ImageSize);
      ardfMRI_denoised_unringed_sdc_b_0(:,:,:,:) = ardfMRI_denoised_unringed_sdc_paired(:,:,1:No_slices,:);
      ardfMRI_denoised_unringed_sdc_b_1(:,:,:,:) = ardfMRI_denoised_unringed_sdc_paired(:,:,No_slices+1:end,:);
      
      niftiwrite(ardfMRI_denoised_unringed_sdc_b_0, 'ardfMRI_denoised_unringed_sdc_b_0', header_4D, 'compressed', true);
      niftiwrite(ardfMRI_denoised_unringed_sdc_b_1, 'ardfMRI_denoised_unringed_sdc_b_1', header_4D, 'compressed', true);
      
      system('gzip  rdfMRI_denoised_unringed_sdc_paired.nii');
      system('gzip ardfMRI_denoised_unringed_sdc_paired.nii');
      
      clear rdfMRI_denoised_unringed_sdc_paired ardfMRI_denoised_unringed_sdc_b_0 ardfMRI_denoised_unringed_sdc_b_1 ardfMRI_denoised_unringed_sdc_paired
      
   end
   
   clear rdfMRI_denoised_unringed_sdc ardfMRI_denoised_unringed_sdc
   
   
   if make_adc == 0
      
      %-removing negative and non-brain values
      %-without replacing NaNs with zeros, thresholding doesnt work
      %-discussed here:
      %-https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1804&L=FSL&D=0&P=232986
      system('fslmaths ardfMRI_denoised_unringed_sdc -nan   ardfMRI_denoised_unringed_sdc');
      system('fslmaths ardfMRI_denoised_unringed_sdc -thr 0 ardfMRI_denoised_unringed_sdc');
      
   else
      
      %-making ADC map
      
      system('fslmaths ardfMRI_denoised_unringed_sdc_b_1 -div ardfMRI_denoised_unringed_sdc_b_0 adc');
      system('fslmaths adc -log adc');
      
      %-dividing by the difference in b-values
      system(['fslmaths adc -div -' num2str(b_val_diff) ' adc']);
      
      system(['fslmerge -tr adc adc ' num2str(2*TR)]);
      
      copyfile('adc.nii.gz', 'adc_incl_weird_voxels.nii.gz');
      
      %-removing negative, very high (>4, arbitrary threshold) and non-brain values
      %-without replacing NaNs with zeros, thresholding doesn't work
      %-discussed here:
      %-https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1804&L=FSL&D=0&P=232986
      system('fslmaths adc -nan adc');
      system('fslmaths adc -thr 0 -uthr 4 adc');
      
      system('fslmaths adc -Tmean adc_mean');
      
      %-the b = 0pt2 volumes can be treated as a proxy for non-DW SE BOLD
      if contains(image_name, '0pt2_1')
         system('fslmaths ardfMRI_denoised_unringed_sdc_b_0 -mas adc_mean map_0pt2');
      end
      
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-registration to MNI space
   %-the later MELODIC and FEAT analyses are in subject's space!
   
   if isempty(dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain1Warp.nii.gz'])))
      
      %-brain extraction from the subject's anatomical image
      %-https://dpaniukov.github.io/2016/06/06/brain-extraction-with-ants.html
      %-https://figshare.com/articles/ANTs_ANTsR_Brain_Templates/915436
      system(['antsBrainExtraction.sh -d 3 ' ...
         ' -a ' fullfile(path_manage, 'scans', [anat_name '.nii.gz']) ...
         ' -e ' fullfile(path_manage_software, 'registration_to_MNI', 'MICCAI2012-Multi-Atlas-Challenge-Data', 'T_template0.nii.gz') ...
         ' -m ' fullfile(path_manage_software, 'registration_to_MNI', 'MICCAI2012-Multi-Atlas-Challenge-Data', 'T_template0_BrainCerebellumProbabilityMask.nii.gz') ...
         ' -f ' fullfile(path_manage_software, 'registration_to_MNI', 'MICCAI2012-Multi-Atlas-Challenge-Data', 'T_template0_BrainCerebellumRegistrationMask.nii.gz') ...
         ' -o ' anat_name '_']);
      
      rmdir([anat_name '_']);
      
      %-registration of the subject's anatomical image to MNI space
      system(['antsRegistrationSyN.sh -d 3 ' ...
         ' -f ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_1mm_brain.nii') ...
         ' -m ' anat_name '_BrainExtractionBrain.nii.gz -o ' anat_name '_anat2mni_brain -n 24']);
      
      %-moving brain extraction and MNI registration files to a folder dedicated to these files for all subjects
      %-important, as I often repeat analyses, and registrations to MNI space take much time
      files_to_copy = dir([anat_name '*']);
      %-the destination folder needs to exist
      mkdir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9)));
      for file_id   = 1:length(files_to_copy)
         movefile(files_to_copy(file_id).name, fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), files_to_copy(file_id).name));
      end
      
   end
   
   %-calculation of the mean of the functional image
   %-this mean map will be registered to subject's anatomical image in next step
   if make_adc == 0
      func_name = 'ardfMRI_denoised_unringed_sdc';
      system('fslmaths ardfMRI_denoised_unringed_sdc.nii.gz     -Tmean func_mean');
   else
      func_name = 'adc';
      system('fslmaths ardfMRI_denoised_unringed_sdc_b_0.nii.gz -Tmean func_mean');
   end
   
   %-https://sourceforge.net/p/advants/discussion/840261/thread/7d71cb7d64/
   %-registration of the subject's functional image to his/her anatomical image
   
   system(['antsRegistration --dimensionality 3 --float 0 ' ...
      ' --output [func2anat, func2anat.nii] ' ...
      ' --interpolation Linear ' ...
      ' --winsorize-image-intensities [0.01,0.99] ' ...
      ' --use-histogram-matching 0 ' ...
      ' --transform Rigid[0.0001] ' ...
      ' --metric MI[' fullfile(path_manage, 'scans', [anat_name '.nii.gz']) ',func_mean.nii.gz,1,32,Regular,0.25]' ...
      ' --convergence [10000x10000x10000x10000,1e-6,10] ' ...
      ' --shrink-factors 4x4x2x1 ' ...
      ' --smoothing-sigmas 2x2x1x0vox']);
   
   %-registration of the subject's functional image to MNI space
   %-combining the two transformations obtained above
   %-only one registration from anatomical image to MNI, per subject
   
   anat2mni_affine_dir   = dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain0GenericAffine.mat']));
   anat2mni_warp_dir     = dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain1Warp.nii.gz']));
   anat2mni_inv_warp_dir = dir(fullfile(path_manage, 'registrations_SyN_subj_anat_to_MNI', image_name(1:9), [anat_name '_anat2mni_brain1InverseWarp.nii.gz']));
   
   if length(anat2mni_affine_dir)~=1 || length(anat2mni_warp_dir)~=1 || length(anat2mni_inv_warp_dir)~=1
      disp(['problem with ' num2str(analysis_id)]);
      continue
   end
   anat2mni_affine   = fullfile(anat2mni_affine_dir(1).folder,   anat2mni_affine_dir(1).name);
   anat2mni_warp     = fullfile(anat2mni_warp_dir(1).folder,     anat2mni_warp_dir(1).name);
   anat2mni_inv_warp = fullfile(anat2mni_inv_warp_dir(1).folder, anat2mni_inv_warp_dir(1).name);
   
   %-bringing the functional image to MNI space
   system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
      ' -i ' func_name '.nii.gz -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v -o func2mni.nii.gz']);
   
   %-bringing maps to subject's space
   
   system(['antsApplyTransforms -r func_mean.nii.gz -n NearestNeighbor ' ...
      ' -i ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
      ' -e 0 -t [func2anat0GenericAffine.mat, 1] [' anat2mni_affine ', 1] ' anat2mni_inv_warp ' -v -o MNI152_T1_2mm_brain_subject_space.nii.gz']);
   
   system(['antsApplyTransforms -r func_mean.nii.gz -n NearestNeighbor ' ...
      ' -i ' atlas_path ...
      ' -e 0 -t [func2anat0GenericAffine.mat, 1] [' anat2mni_affine ', 1] ' anat2mni_inv_warp ' -v -o NMM_atlas_subject_space.nii.gz']);
   
   system(['antsApplyTransforms -r func_mean.nii.gz -n NearestNeighbor ' ...
      ' -i ' atlas_w_JHU_path ...
      ' -e 0 -t [func2anat0GenericAffine.mat, 1] [' anat2mni_affine ', 1] ' anat2mni_inv_warp ' -v -o NMM_JHU_atlas_subject_space.nii.gz']);
   
   for VOI = {'all', 'visual', 'motor', 'somatosensory'}
      mask_path = fullfile(path_manage_software, 'brain_parcellation', ['mask_VOI_' VOI{1} '.nii']);
      system(['antsApplyTransforms -r func_mean.nii.gz -n NearestNeighbor ' ...
         ' -i ' mask_path ' -e 0 -t [func2anat0GenericAffine.mat, 1] [' anat2mni_affine ', 1] ' anat2mni_inv_warp ' -v -o mask_VOI_' VOI{1} '_subject_space.nii.gz']);
   end
   
   %-increasing/dilating the masks, to account for registration imperfections
   for VOI = {'all', 'visual', 'motor', 'somatosensory'}
      system(['fslmaths mask_VOI_' VOI{1} '_subject_space         -bin                    mask_VOI_' VOI{1} '_subject_space_dilated']);
      system(['fslmaths mask_VOI_' VOI{1} '_subject_space_dilated -kernel sphere 3 -fmean mask_VOI_' VOI{1} '_subject_space_dilated']);
      system(['fslmaths mask_VOI_' VOI{1} '_subject_space_dilated -bin                    mask_VOI_' VOI{1} '_subject_space_dilated']);
   end
   
   %-making brain mask
   copyfile('mask_VOI_all_subject_space_dilated.nii.gz', 'brain_mask.nii.gz');
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-preparations for the resting state and task analyses, which follow
   
   if make_adc == 0
      replace_TR           = TR;
      replace_No_volumes   = num2str(No_volumes);
      replace_total_voxels = num2str(prod(header_4D.ImageSize));
      replace_data_part    = 'ardfMRI_denoised_unringed_sdc';
   else
      replace_TR           = 2*TR;
      replace_No_volumes   = num2str(floor(No_volumes/2));
      replace_total_voxels = num2str(prod(header_4D.ImageSize));
      replace_data_part    = 'adc';
   end
   
   replace_data            = fullfile(path_manage, 'analyses_output', image_name, replace_data_part);
   replace_mask            = fullfile(path_manage, 'analyses_output', image_name, 'brain_mask.nii.gz');
   
   %-the above mask didnt help!
   %-https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=ind02&L=FSL&P=R53261https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=ind02&L=FSL&P=R53261
   %-Brain/background threshold needs to be set to 0!
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-analysis of resting state data (MELODIC)
   
   %-copy the MELODIC's design file
   copyfile(fullfile(path_manage_software, 'design_rest.fsf'), fullfile(path_manage, 'analyses_output', image_name, 'design_rest.fsf'));
   
   %-':' used in sed instead of '\', due to '\' being used in paths
   system(['sed -i ''s:replace_TR:'             num2str(replace_TR)    ':g'' design_rest.fsf']);
   system(['sed -i ''s:replace_No_volumes:'     replace_No_volumes     ':g'' design_rest.fsf']);
   system(['sed -i ''s:replace_total_voxels:'   replace_total_voxels   ':g'' design_rest.fsf']);
   system(['sed -i ''s:replace_data:'           replace_data           ':g'' design_rest.fsf']);
   system(['sed -i ''s:replace_dim_estimation:' replace_dim_estimation ':g'' design_rest.fsf']);
   system(['sed -i ''s:replace_output_comps:'   replace_output_comps   ':g'' design_rest.fsf']);
   system(['sed -i ''s:replace_mask:'           replace_mask           ':g'' design_rest.fsf']);
   
   %-MELODIC is some kind of FEAT
   system('feat design_rest.fsf');
   
   %-bringing the filtered image to MNI space (for the resting state analysis)
   system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
      ' -i ' fullfile([replace_data '.ica'], 'filtered_func_data.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
      ' -o ' fullfile([replace_data '.ica'], 'filtered_func_data2mni.nii')]);
   
   %-the b = 0pt2 volumes can be treated as a proxy for non-DW SE BOLD
   if contains(image_name, '0pt2_1')
      
      %-adjust the MELODIC's design file
      %-':' used in sed instead of '\', due to '\' being used in paths
      system('sed -i ''s:/adc:/map_0pt2:g'' design_rest.fsf');
      
      %-MELODIC is some kind of FEAT
      system('feat design_rest.fsf');
      
      %-bringing the filtered image to MNI space (for the resting state analysis)
      system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
         ' -i ' fullfile('map_0pt2.ica', 'filtered_func_data.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
         ' -o ' fullfile('map_0pt2.ica', 'filtered_func_data2mni.nii')]);
      
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-analysis of task data (FEAT)
   
   if task_analysis == 1
      
      %-copy the FEAT's design file
      copyfile(fullfile(path_manage_software, 'design_task.fsf'), fullfile(path_manage, 'analyses_output', image_name, 'design_task.fsf'));
      pause(2);
      
      %-creating a '.txt' file with the experimental design
      
      %-m-sequence
      if contains(image_name, '_m_seq_')
         
         stim_onsets  = readtable(fullfile(path_manage_software, 'mseq_stim_onsets.txt'), 'ReadVariableNames', false);
         stim_onsets  = table2array(stim_onsets);
         exper_design = zeros(128, 3);
         exper_design(:,1) = stim_onsets;
         exper_design(:,3) = 1;
         
      else
         
         exper_design = zeros(1000, 3);
         
         for i = 1:200
            exper_design(i, 1) = (i-1)*(rest_block+task_block)-No_dummies_out*TR+rest_block;
            %-duration 0, due to FIR modelling
            exper_design(i, 2) = task_block;
            exper_design(i, 3) = 1;
            if exper_design(i, 1) > str2double(replace_No_volumes)*replace_TR
               exper_design = exper_design(1:i-1,:);
               break
            end
         end
         
      end
      
      writetable(table(exper_design), 'exper_design.txt', 'Delimiter', ' ', 'WriteVariableNames', false);
      
      
      %-':' used in sed instead of '\', due to '\' being used in paths
      system(['sed -i ''s:replace_TR:'                num2str(replace_TR)               ':g'' design_task.fsf']);
      system(['sed -i ''s:replace_No_volumes:'        replace_No_volumes                ':g'' design_task.fsf']);
      system(['sed -i ''s:replace_total_voxels:'      replace_total_voxels              ':g'' design_task.fsf']);
      system(['sed -i ''s:replace_data:'              replace_data                      ':g'' design_task.fsf']);
      system(['sed -i ''s:replace_design:'            fullfile(pwd, 'exper_design.txt') ':g'' design_task.fsf']);
      system(['sed -i ''s:replace_FIR_window_length:' num2str(14*2*TR)                  ':g'' design_task.fsf']);
      
      system('feat design_task.fsf');
      
      feat_folder = dir('*.feat');
      
      %-for m-sequence the last parameters do not make sense!
      plot_diagnostic_plots(fullfile(feat_folder.folder, feat_folder.name, 'stats'), 'res4d', 'brain_mask', replace_TR, 1/100, 1/(rest_block+task_block), 1/(rest_block+task_block));
      
      %-bringing the filtered image to MNI space (for the task analysis)
      system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
         ' -i ' fullfile([replace_data '.feat'], 'filtered_func_data.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
         ' -o ' fullfile([replace_data '.feat'], 'filtered_func_data2mni.nii')]);
      
      %-bringing zfstat1 to MNI space
      system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
         ' -i ' fullfile([replace_data '.feat'], 'stats', 'zfstat1.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
         ' -o ' fullfile([replace_data '.feat'], 'stats', 'zfstat1_MNI.nii.gz')]);
      
      %-the b = 0pt2 volumes can be treated as a proxy for non-DW SE BOLD
      if contains(image_name, '0pt2_1')
         
         %-adjust the FEAT's design file
         %-':' used in sed instead of '\', due to '\' being used in paths
         system('sed -i ''s:/adc:/map_0pt2:g'' design_task.fsf');
         
         system('feat design_task.fsf');
         
         %-bringing the filtered image to MNI space (for the task analysis)
         system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
            ' -i ' fullfile('map_0pt2.feat', 'filtered_func_data.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
            ' -o ' fullfile('map_0pt2.feat', 'filtered_func_data2mni.nii')]);
         
         %-bringing zfstat1 to MNI space
         system(['antsApplyTransforms -r ' fullfile(path_manage_software, 'registration_to_MNI', 'MNI152_T1_2mm_brain.nii') ...
            ' -i ' fullfile('map_0pt2.feat', 'stats', 'zfstat1.nii.gz') ' -e 3 -t ' anat2mni_warp ' ' anat2mni_affine ' func2anat0GenericAffine.mat -v ' ...
            ' -o ' fullfile('map_0pt2.feat', 'stats', 'zfstat1_MNI.nii.gz')]);
         
      end
      
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-summarising SNRs, done after calculating SNR maps as now we have a brain mask
   
   brain_mask               = niftiread('brain_mask.nii.gz');
   SNR_before_PCA_denoising = niftiread('SNR/SNR_before_PCA_denoising.nii.gz');
   SNR_after_PCA_denoising  = niftiread('SNR/SNR_after_PCA_denoising.nii.gz');
   
   fprintf(fileID, ['adc: ' num2str(make_adc) '\n\n']);
   %-only for three images there are interleaved b=0,1 values
   %-for the other three images the means should be very similar
   SNR_before_PCA_denoising_b_0 = SNR_before_PCA_denoising(:,:,:,1);
   SNR_before_PCA_denoising_b_1 = SNR_before_PCA_denoising(:,:,:,2);
   fprintf(fileID, ['before PCA: 1st volume: ' num2str(mean(SNR_before_PCA_denoising_b_0(brain_mask==1))) '\n']);
   fprintf(fileID, ['before PCA: 2nd volume: ' num2str(mean(SNR_before_PCA_denoising_b_1(brain_mask==1))) '\n']);
   SNR_after_PCA_denoising_b_0  = SNR_after_PCA_denoising(:,:,:,1);
   SNR_after_PCA_denoising_b_1  = SNR_after_PCA_denoising(:,:,:,2);
   fprintf(fileID, ['after PCA: 1st volume: '  num2str(mean(SNR_after_PCA_denoising_b_0(brain_mask==1)))  '\n']);
   fprintf(fileID, ['after PCA: 2nd volume: '  num2str(mean(SNR_after_PCA_denoising_b_1(brain_mask==1)))  '\n']);
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-assessing image smoothness
   
   fprintf(fileID, '\n dfMRI_raw:                    \n');
   system('smoothest -z dfMRI_raw                    -m brain_mask >> log.txt');
   fprintf(fileID, '\n dfMRI_denoised:               \n');
   system('smoothest -z dfMRI_denoised               -m brain_mask >> log.txt');
   fprintf(fileID, '\n dfMRI_denoised_unringed:      \n');
   system('smoothest -z dfMRI_denoised_unringed      -m brain_mask >> log.txt');
   fprintf(fileID, '\n dfMRI_denoised_unringed_sdc:  \n');
   system('smoothest -z dfMRI_denoised_unringed_sdc  -m brain_mask >> log.txt');
   fprintf(fileID, '\n rdfMRI_denoised_unringed_sdc: \n');
   system('smoothest -z rdfMRI_denoised_unringed_sdc -m brain_mask >> log.txt');
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-making diagnostic plots
   
   plot_diagnostic_plots(pwd,                    'dfMRI_raw',               'brain_mask', TR, 0.01);
   plot_diagnostic_plots(pwd,                    'dfMRI_denoised',          'brain_mask', TR, 0.01);
   plot_diagnostic_plots(pwd,                    'PCA_denoising_residuals', 'brain_mask', TR, 0.01);
   plot_diagnostic_plots(fullfile(pwd, 'SNR'),   'dfMRI_raw_hpf',           'brain_mask', TR, 0.01);
   plot_diagnostic_plots(fullfile(pwd, 'SNR'),   'dfMRI_denoised_hpf',      'brain_mask', TR, 0.01);
   plot_diagnostic_plots([replace_data '.ica'],  'filtered_func_data2mni',  'none',       TR, 0.01);
   
   system(['gunzip ' fullfile([replace_data '.ica'],  'filtered_func_data2mni.nii.gz')]);
   
   if task_analysis == 1
      plot_diagnostic_plots([replace_data '.feat'], 'filtered_func_data2mni', 'none', TR, 0.01);
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %-plotting correlations between the ROI mean time series
   
   for map_0pt2_yn = 0:1
      
      if map_0pt2_yn == 1 && ~contains(image_name, '0pt2_1')
         continue
      end
      
      if map_0pt2_yn == 0
         fig_name_suffix_2 = '';
      else
         fig_name_suffix_2 = '_map_0pt2';
      end
      
      figure('visible', 'off');
      
      atlas_parcellation = niftiread('NMM_atlas_subject_space.nii.gz');
      
      for clean_yes = 0:1
         
         if task_analysis == 0
            if clean_yes == 0
               if map_0pt2_yn == 0
                  func_map     = niftiread(fullfile([replace_data '.ica'], 'filtered_func_data.nii.gz'));
               else
                  func_map     = niftiread(fullfile('map_0pt2.ica',        'filtered_func_data.nii.gz'));
               end
               fig_name_suffix = '';
            else
               if     map_0pt2_yn == 0 && contains(image_name, '_diff_')
                  ICA_suffix = 'diff_rs';
               elseif map_0pt2_yn == 0 && contains(image_name, '_se_')
                  ICA_suffix = 'se_rs';
               else
                  ICA_suffix = 'se_proxy';
               end
               cleaned_map  = dir(fullfile(path_manage, 'ICA', ICA_suffix, image_name, '*ic40_clean.nii.gz'));
               if length(cleaned_map) ~= 1
                  disp(['problems with cleaned image for ' image_name]);
                  continue
               end
               func_map        = niftiread(fullfile(cleaned_map.folder, cleaned_map.name));
               fig_name_suffix = '_clean';
            end
         else
            if clean_yes == 0
               if map_0pt2_yn == 0
                  func_map     = niftiread(fullfile([replace_data '.feat'], 'stats', 'res4d.nii.gz'));
               else
                  func_map     = niftiread(fullfile('map_0pt2.feat',        'stats', 'res4d.nii.gz'));
               end
               fig_name_suffix = '';
            else
               %-manual ICA cleaning was done only for RS datasets
               continue
            end
         end
         
         ROI_means  = zeros(size(func_map, 4), length(ROIs));
         
         for ROI_id = 1:length(ROIs)
            voxels_ids        = find(atlas_parcellation==ROIs(ROI_id));
            for time_id       = 1:size(func_map, 4)
               func_map_1_vol = func_map(:,:,:,time_id);
               %-only considering the ROI if > 10 voxels available (arbitrary threshold)
               if length(voxels_ids) > 10
                  ROI_means(time_id, ROI_id) = mean(func_map_1_vol(voxels_ids));
               end
            end
         end
         
         save(['ROI_means' fig_name_suffix fig_name_suffix_2 '.mat'], 'ROI_means');
         
         %-performing global signal regression (GSR)
         
         global_signal = zeros(1, size(func_map, 4));
         
         for time_id       = 1:size(func_map, 4)
            func_map_1_vol = func_map(:,:,:,time_id);
            global_signal(time_id) = mean(func_map_1_vol(atlas_parcellation>0));
         end
         
         for ROI_id = 1:length(ROIs)
            %-for older versions of matlab (<2020) .Residuals doesnt work in the same line as fitlm
            %-sometimes in ROI_means there are NaNs, then fitlm complains about rank deficiency
            gs_lm = fitlm(global_signal, ROI_means(:, ROI_id));
            ROI_means(:, ROI_id) = table2array(gs_lm.Residuals(:,1));
         end
         
         save(['ROI_means' fig_name_suffix fig_name_suffix_2 '_GSR.mat'], 'ROI_means');
         
         subplot(3, 3, 3+clean_yes*3+1);
         non_NA_indices = find((ROI_means(1,:)~=0) & (~isnan(ROI_means(1,:))));
         imagesc(corr(ROI_means(:, non_NA_indices)));
         caxis([-1 1]);
         colormap('jet');
         ax            = gca;
         ax.YTick      = 1:length(non_NA_indices);
         ax.XTick      = 1:length(non_NA_indices);
         %-the spaces in YTickLabels appear when exporting to pdf, it is a known problem in matlab
         %-I can get rid of it changing the interpreter, but then I cant use the Helvetica font
         %-'TickLabelInterpreter', 'latex', 'FontName', 'Helvetica'
         ax.YTickLabel = ROIs_t(non_NA_indices);
         ax.XTickLabel = ROIs_t(non_NA_indices);
         xtickangle(90);
         set(ax, 'xaxisLocation', 'top');
         set(ax, 'FontSize', 0.5);
         
         ROI_means_corrs = reshape(corr(ROI_means), 1, []);
         ROI_means_corrs = abs(ROI_means_corrs);
         [M, I]          = sort(ROI_means_corrs, 'descend');
         ROIs_with_high_corrs = zeros(1,28);
         for corr_highest_id = 1:14
            id_to_decode = I(max(find(M==1)) + corr_highest_id);
            column_id    = fix((id_to_decode-1)/length(ROIs)) + 1;
            row_id       = id_to_decode - (column_id-1)*length(ROIs);
            ROIs_with_high_corrs(2*corr_highest_id-1) = row_id;
            ROIs_with_high_corrs(2*corr_highest_id)   = column_id;
         end
         
         ROIs_with_high_corrs = unique(ROIs_with_high_corrs);
         
         subplot(3, 3, 3+clean_yes*3+2);
         imagesc(corr(ROI_means(:, ROIs_with_high_corrs)));
         caxis([-1 1]);
         colormap('jet');
         ax            = gca;
         ax.YTick      = 1:length(ROIs_with_high_corrs);
         ax.XTick      = 1:length(ROIs_with_high_corrs);
         ax.YTickLabel = ROIs_t(ROIs_with_high_corrs);
         ax.XTickLabel = ROIs_t(ROIs_with_high_corrs);
         xtickangle(90);
         set(ax, 'xaxisLocation', 'top');
         set(gca, 'FontSize', 2);
         if clean_yes == 0
            title([fig_title_text ':   not cleaned'], 'FontSize', 6);
         else
            title('cleaned', 'FontSize', 6);
         end
         
         %{
	from the Neuromorphometrics atlas:
	23: Right Accumbens Area,
	30: Left Accumbens Area,
	59: Right Thalamus Proper,
	60: Left Thalamus Proper,
	100: Right ACgG anterior cingulate gyrus,
	101: Left ACgG anterior cingulate gyrus,
	166: Right PCgG posterior cingulate gyrus,
	167: Left PCgG posterior cingulate gyrus,
	108: Right Calc calcarine cortex,
	109: Left Calc calcarine cortex,
	182: Right PrG precentral gyrus,
	183: Left PrG precentral gyrus,
	192: Right SMC supplementary motor cortex,
	193: Left SMC supplementary motor cortex
         %}
         
         ROIs_for_imagesc = [23 30 59 60 100 101 166 167 108 109 182 183 192 193];
         ROIs_list = zeros(1, length(ROIs_for_imagesc));
         for i = 1:length(ROIs_for_imagesc)
            ROIs_list(i) = find(ROIs_for_imagesc(i)==ROIs);
         end
         
         subplot(3, 3, 3+clean_yes*3+3);
         imagesc(corr(ROI_means(:, ROIs_list)));
         caxis([-1 1]);
         colormap('jet');
         ax            = gca;
         ax.YTick      = 1:length(ROIs_list);
         ax.XTick      = 1:length(ROIs_list);
         ax.YTickLabel = ROIs_t(ROIs_list);
         ax.XTickLabel = ROIs_t(ROIs_list);
         xtickangle(90);
         set(ax, 'xaxisLocation', 'top');
         cb = colorbar;
         cb.Position(1) = 0.93; % or + [left, bottom, width, height] to place it where you want
         set(gca, 'FontSize', 2);
         
      end
      
      h = gcf;
      set(h, 'PaperOrientation', 'landscape');
      set(h, 'Units', 'Inches');
      pos = get(h, 'Position');
      set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
      if analysis_id < 10
         zeros_str = '0';
      else
         zeros_str = '';
      end
      print(['ROI_correlations_adc_' num2str(make_adc) '_' zeros_str num2str(analysis_id) fig_name_suffix_2], '-dpdf', '-r0');
      system(['pdfcrop ROI_correlations_adc_' num2str(make_adc) '_' zeros_str num2str(analysis_id) fig_name_suffix_2 '.pdf ROI_correlations_adc_' num2str(make_adc) '_' zeros_str num2str(analysis_id) fig_name_suffix_2 '.pdf']);
      
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if task_analysis == 1
      
      %-if the entire pipeline is run, 'exper_design' is in the memory
      exper_design = readtable('exper_design.txt', 'Delimiter', ' ', 'ReadVariableNames', false);
      exper_design = table2array(exper_design);
      
      %-for the 1st image, F test contrast error, TR_ADC = 2s (not 2.07!), design 12+18 (reversed than the majority)
      %-for the 2nd image, too few volumes for FIR
      %-these were only pilots
      if strcmp(image_name, 'pilot_EF_Dfmri_64cx_sms2_2iso_diff_600dir_b0_1_task_12s_off_18s_on_20200312162637_9') || ...
            strcmp(image_name, 'pilot_WO2_cmrr_mbep2d_diff_b0_1_task_18s_off_2s_on_20200520100309_7')
         continue
      end
      
      system(' gunzip func2mni.nii.gz');
      system(['gunzip ' replace_data_part '.nii.gz']);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %-running task GLM through SPM, both in subject's space and in MNI space
      %-MELODIC and FEAT analyses (above) are in subject's space
      
      %-this is awkward, 'replace_TR' is sometimes 2*1.035...   !!!!!!!
      if replace_TR == 2
         No_FIR_bins_SPM = 15;
      else
         No_FIR_bins_SPM = 30;
      end
      
      for space = {'subject', 'MNI'}
         
         scans  = strings(str2double(replace_No_volumes), 1);
         
         for volume_id  = 1:str2double(replace_No_volumes)
            
            if strcmp(space{1}, 'MNI')
               scans{volume_id} = [fullfile(pwd, 'func2mni.nii')             ',' num2str(volume_id)];
            else
               scans{volume_id} = [fullfile(pwd, [replace_data_part '.nii']) ',' num2str(volume_id)];
            end
            
         end
         
         clear jobs;
         
         %-model specification and estimation
         jobs{1}.stats{1}.fmri_spec.dir          = cellstr(fullfile(pwd, ['SPM_task_GLM_' space{1} '_space']));
         jobs{1}.stats{1}.fmri_spec.timing.units = 'secs';
         jobs{1}.stats{1}.fmri_spec.timing.RT    = replace_TR;
         
         %-with standard 16/8 problems with the contrast for the ADC F-test
         jobs{1}.stats{1}.fmri_spec.timing.fmri_t  = 32;
         jobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = 17;
         jobs{1}.stats{1}.fmri_spec.sess.scans   = cellstr(scans);
         
         %-specifying the experimental design
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).name     = 'checkerboard';
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).onset    = exper_design(:,1);
         jobs{1}.stats{1}.fmri_spec.sess.cond(1).duration = 0;
         
         jobs{1}.stats{1}.fmri_spec.bases.fir.length      = 30;
         jobs{1}.stats{1}.fmri_spec.bases.fir.order       = No_FIR_bins_SPM;
         
         %-no brain mask calculated
         jobs{1}.stats{1}.fmri_spec.mthresh   = -Inf;
         jobs{1}.stats{1}.fmri_spec.mask      = {''};
         
         jobs{1}.stats{2}.fmri_est.spmmat(1)  = cellstr(['SPM_task_GLM_' space{1} '_space/SPM.mat']);
         
         %-save residuals
         jobs{1}.stats{2}.fmri_est.write_residuals = 1;
         
         jobs{1}.stats{3}.con.spmmat          = cellstr(['SPM_task_GLM_' space{1} '_space/SPM.mat']);
         jobs{1}.stats{3}.con.consess{1}.fcon = struct('name', 'activation', 'convec', eye(No_FIR_bins_SPM), 'sessrep', 'none');
         
         spm_jobman('run', jobs);
         
         %-merging SPM's GLM residuals to one file
         cd(['SPM_task_GLM_' space{1} '_space']);
         res_all    = dir('Res_*.nii');
         res_all    = {res_all.name};
         res_all_4d = char(res_all);
         %-TR can not be specified...
         
         spm_file_merge(res_all_4d, 'res4d.nii');
         for res_id = 1:str2double(replace_No_volumes)
            delete(char(res_all(res_id)));
         end
         
         cd ..
         
         plot_diagnostic_plots(fullfile(pwd, ['SPM_task_GLM_' space{1} '_space']), 'res4d', 'none', TR, 1/128);
         
         %-transforming the F-statistic map to a z-statistic map
         load(['SPM_task_GLM_' space{1} '_space/SPM.mat']);
         system(['ftoz -zout SPM_task_GLM_' space{1} '_space/zfstat1 SPM_task_GLM_' space{1} '_space/spmF_0001.nii ' num2str(No_FIR_bins_SPM+1) ' ' num2str(SPM.xX.erdf)]);
         clear SPM
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %-calculating the mean time-courses and response functions
      
      %-calculating the VOI/PCA time-courses
      %-currently, the mean time-courses and response functions are based on the mean signal;
      %-however, because there might be both positive and negative components in the signal,
      %-using 1st PCA component instead of the mean would be better;
      %-I am not sure whether this calculation of the VOI/PCA time-course is correct!
      
      for VOI     = {'visual', 'motor', 'somatosensory'}
         
         %-in 'extra_functions' there is Rik's version of 'spm_regions'
         %-important for avoiding extreme values in 1st and last volumes
         %-adjust '0' also led to extreme values in 1st and last volumes
         
         mask_path = fullfile(path_manage_software, 'brain_parcellation', ['mask_VOI_' VOI{1} '.nii']);
         clear jobs;
         jobs{1}.util{1}.voi.spmmat                = cellstr(fullfile(pwd, 'SPM_task_GLM_subject_space', 'SPM.mat'));
         jobs{1}.util{1}.voi.adjust                = NaN;
         jobs{1}.util{1}.voi.session               = 1;
         jobs{1}.util{1}.voi.name                  = VOI{1};
         jobs{1}.util{1}.voi.roi{1}.mask.image     = cellstr(mask_path);
         jobs{1}.util{1}.voi.roi{1}.mask.threshold = 0.5;
         jobs{1}.util{1}.voi.expression            = 'i1';
         spm_jobman('run', jobs);
         
         %-the VOI calculation seems to open a window
         close all;
         
      end
      
      for map_0pt2_yn = 0:1
         
         if map_0pt2_yn == 1 && ~contains(image_name, '0pt2_1')
            continue
         end
         
         if map_0pt2_yn == 0
            fig_name_suffix_2 = '';
         else
            fig_name_suffix_2 = '_map_0pt2';
         end
         
         figure('visible', 'off', 'pos', [0 0 800 800]);
         
         set(0, 'DefaultAxesTitleFontWeight', 'normal');
         
         if make_adc == 0
            func_map = niftiread(fullfile([replace_data '.feat'], 'filtered_func_data.nii.gz'));
         elseif map_0pt2_yn == 0
            if exist('adc.nii', 'file') == 2
               func_map = niftiread('adc.nii');
            else
               func_map = niftiread('adc.nii.gz');
            end
         else
            func_map = niftiread(fullfile('map_0pt2.feat',        'filtered_func_data.nii.gz'));
         end
         
         if map_0pt2_yn == 0
            cd(fullfile([replace_data '.feat'], 'stats'));
         else
            cd(fullfile('map_0pt2.feat',        'stats'));
         end
         
         %-arbitrary thresholds
         system('cluster -i zfstat1 -t 3.1 --osize=cluster_size > cluster_info.txt');
         
         cluster_size = niftiread('cluster_size.nii.gz');
         
         cd ../..
         
         if make_adc == 1
            system('fslmaths adc -sub adc_incl_weird_voxels adc_weird_voxels');
            system('fslmaths adc_weird_voxels -abs adc_weird_voxels');
            system('fslmaths adc_weird_voxels -bin adc_weird_voxels');
            system('fslmaths adc_weird_voxels -Tmean adc_weird_voxels_perc');
            adc_weird_voxels_perc = niftiread('adc_weird_voxels_perc.nii.gz');
         end
         
         for VOI     = {'all', 'visual', 'motor', 'somatosensory'}
            
            mask_VOI = niftiread(['mask_VOI_' VOI{1} '_subject_space_dilated.nii.gz']);
            
            %-for 'clustered', I first used 10, then 4
            mask_VOI(cluster_size < 1) = 0;
            
            if make_adc == 1
               mask_VOI(adc_weird_voxels_perc > 0.05) = 0;
            end
            
            %-calculating the mean time-course
            
            voxels_ids                   = find(mask_VOI > 0);
            time_course_mean             = zeros(1, size(func_map, 4));
            for time_id                  = 1:size(func_map, 4)
               func_map_1_vol            = func_map(:,:,:,time_id);
               time_course_mean(time_id) = mean(func_map_1_vol(voxels_ids));
            end
            
            % save(['time_course_mean_' VOI{1} fig_name_suffix_2 '_clustered'], 'time_course_mean');
            save(['time_course_mean_' VOI{1} fig_name_suffix_2], 'time_course_mean');
            
            %-plotting the time-courses
            
            ts = time_course_mean;
            
            if strcmp(VOI{1}, 'all')
               subplot(3, 4, 1);
               title_text = {['all brain: '        num2str(sum(func_map_1_vol(voxels_ids)>0)) ' vox'], ' ', ' '};
            elseif strcmp(VOI{1}, 'visual')
               subplot(3, 4, 2);
               title_text = {['visual cortex: '    num2str(sum(func_map_1_vol(voxels_ids)>0)) ' vox'], ' ', ' '};
            elseif strcmp(VOI{1}, 'motor')
               subplot(3, 4, 3);
               title_text = {['motor cortex: '     num2str(sum(func_map_1_vol(voxels_ids)>0)) ' vox'], ' ', ' '};
            else
               subplot(3, 4, 4);
               title_text = {['somatosensory c: '  num2str(sum(func_map_1_vol(voxels_ids)>0)) ' vox'], ' ', ' '};
            end
            
            plot((1:str2double(replace_No_volumes))*replace_TR, ts);
            for i = 1:size(exper_design, 1)
               patch(exper_design(i,1)+[0 task_block task_block 0], [min(ylim) min(ylim) max(ylim) max(ylim)], [200 200 200]/255, 'EdgeColor', 'none');
            end
            hold on;
            plot((1:str2double(replace_No_volumes))*replace_TR, ts);
            if sum(ts) ~= 0 && sum(isnan(ts)) == 0
               ylim([min(ts) max(ts)]);
            end
            htitle = title(title_text, 'FontSize', 9);
            xlabel('time [s]');
            xlim([0 str2double(replace_No_volumes)*replace_TR]);
            
            if sum(isnan(time_course_mean)) == 0 && sum(time_course_mean) ~= 0
               
               %-calculating and plotting the response functions: dRF and HRF
               
               outliers_low  = find(time_course_mean < (mean(time_course_mean) - 2*std(time_course_mean)));
               outliers_high = find(time_course_mean > (mean(time_course_mean) + 2*std(time_course_mean)));
               outliers      = [outliers_low outliers_high];
               non_outliers  = setdiff(1:str2double(replace_No_volumes), outliers);
               
               load SPM_task_GLM_subject_space/SPM
               
               SPM.xX.X_wo_outliers         = SPM.xX.X        (non_outliers, :);
               time_course_mean_wo_outliers = time_course_mean(non_outliers);
               
               mdl = fitlm(SPM.xX.X_wo_outliers, time_course_mean_wo_outliers, 'Intercept', false);
               
               clear SPM
               
               if strcmp(VOI{1}, 'all')
                  subplot(3, 4, 5);
               elseif strcmp(VOI{1}, 'visual')
                  subplot(3, 4, 6);
               elseif strcmp(VOI{1}, 'motor')
                  subplot(3, 4, 7);
               else
                  subplot(3, 4, 8);
               end
               
               if replace_TR == 2
                  bins_middle_points =   1 : 2 : 29;
               else
                  bins_middle_points = 0.5 : 1 : 29.5;
               end
               
               %-the last coefficient corresponds to the intercept
               response = table2array(mdl.Coefficients(1:end-1, 1));
               bar(bins_middle_points, response);
               patch([0 task_block task_block 0], [min(ylim) min(ylim) max(ylim) max(ylim)], [200 200 200]/255, 'EdgeColor', 'none');
               hold on;
               bar(bins_middle_points, response);
               %-to get a different color than for the time-courses
               bar(bins_middle_points, response);
               xlim([-8 30]);
               ylim([min(0, min(response)) max(response)]);
               xticks([0 10 20 30]);
               xticklabels({'0','10','20','30'});
               
               title({' '; ['volume outliers: ' num2str(length(outliers))]; ' '}, 'FontSize', 9);
               xlabel('time [s]');
               
               %-averaging signal across trials (and plotting), checking if similar response function as with FIR
               
               time_course_mean_wo_outliers = time_course_mean;
               time_course_mean_wo_outliers(outliers) = NaN;
               
               %-interpolation doesnt work for first and last elements
               %-they cannot be NaN
               time_course_mean_wo_outliers(1)   = time_course_mean_wo_outliers(non_outliers(1));
               time_course_mean_wo_outliers(end) = time_course_mean_wo_outliers(non_outliers(end));
               
               %-https://ch.mathworks.com/matlabcentral/answers/34346-interpolating-nan-s
               nanx = isnan(time_course_mean_wo_outliers);
               t    = 1:numel(time_course_mean_wo_outliers);
               time_course_mean_wo_outliers(nanx) = interp1(t(~nanx), time_course_mean_wo_outliers(~nanx), t(nanx));
               
               time          = (0.5    : 1     : (str2double(replace_No_volumes)-0.5))*replace_TR;
               time_accurate = time(1) : 0.001 : time(end);
               
               time_course_interpolated = interp1q(time', time_course_mean_wo_outliers', time_accurate');
               
               response          = zeros(38001, 1);
               No_of_repetitions = zeros(38001, 1);
               
               for trial_id  = 1:length(exper_design(:,1))
                  
                  %-multiplied with 1000 as the desired temporal resolution is 0.001s
                  %-we would like to also see 8s before the stimulus onset
                  time_point_1st  = round(exper_design(trial_id,1)*1000) - 1000*time_accurate(1) - 8000;
                  time_point_last = time_point_1st + 1000*(30+8);
                  
                  %-for the last trials, perhaps not the entire task+rest blocks visible
                  if time_point_last > length(time_course_interpolated)
                     time_point_last = length(time_course_interpolated);
                  end
                  
                  %-for the m-sequence, not possible to get the 8s before the 1st and 2nd stimulus onset
                  if time_point_1st < 1 && (trial_id == 1 || trial_id == 2)
                     continue;
                  end
                  
                  No_time_points    = time_point_last - time_point_1st + 1;
                  
                  %-indices are required to be integers
                  response(int64(1:No_time_points)) = response(int64(1:No_time_points)) + time_course_interpolated(int64(time_point_1st : time_point_last));
                  No_of_repetitions(int64(1:No_time_points)) = No_of_repetitions(int64(1:No_time_points)) + 1;
                  
               end
               
               response = response./No_of_repetitions;
               
               % save(['response_' VOI{1} fig_name_suffix_2 '_clustered'], 'response');
               save(['response_' VOI{1} fig_name_suffix_2], 'response');
               
               if strcmp(VOI{1}, 'all')
                  subplot(3, 4, 9);
               elseif strcmp(VOI{1}, 'visual')
                  subplot(3, 4, 10);
               elseif strcmp(VOI{1}, 'motor')
                  subplot(3, 4, 11);
               else
                  subplot(3, 4, 12);
               end
               
               plot(-8:0.001:30, response);
               patch([0 task_block task_block 0], [min(ylim) min(ylim) max(ylim) max(ylim)], [200 200 200]/255, 'EdgeColor', 'none');
               hold on;
               plot(-8:0.001:30, response);
               xlim([-8 30]);
               ylim([min(response) max(response)]);
               xticks([0 10 20 30]);
               xticklabels({'0','10','20','30'});
               
               xlabel('time [s]');
               
               title({' '; ['volume outliers: ' num2str(length(outliers))]; ' '}, 'FontSize', 9);
               
            end
            
         end
         
         %-adding figure titles
         sgtitle([fig_title_text strrep(fig_name_suffix_2, '_', ' ')], 'FontSize', 10);
         
         annotation('textbox', [0 0.815 1 0.06], 'String', 'time courses',                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
         annotation('textbox', [0 0.565 1 0.06], 'String', 'FIR response functions',      'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
         annotation('textbox', [0 0.284 1 0.06], 'String', 'averaged response functions', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);
         
         if analysis_id < 10
            % print_to_svg_to_pdf(['time_courses_and_response_functions_0' num2str(analysis_id) fig_name_suffix_2 '_clustered']);
            print_to_svg_to_pdf(['time_courses_and_response_functions_0' num2str(analysis_id) fig_name_suffix_2]);
         else
            % print_to_svg_to_pdf(['time_courses_and_response_functions_'  num2str(analysis_id) fig_name_suffix_2 '_clustered']);
            print_to_svg_to_pdf(['time_courses_and_response_functions_'  num2str(analysis_id) fig_name_suffix_2]);
         end
         
      end
      
      plot_diagnostic_plots(pwd, 'time_course_mean_all',           'none', replace_TR, 0.01, 1/(rest_block+task_block), 1/(rest_block+task_block));
      plot_diagnostic_plots(pwd, 'time_course_mean_visual',        'none', replace_TR, 0.01, 1/(rest_block+task_block), 1/(rest_block+task_block));
      plot_diagnostic_plots(pwd, 'time_course_mean_motor',         'none', replace_TR, 0.01, 1/(rest_block+task_block), 1/(rest_block+task_block));
      plot_diagnostic_plots(pwd, 'time_course_mean_somatosensory', 'none', replace_TR, 0.01, 1/(rest_block+task_block), 1/(rest_block+task_block));
      
      %-this will not work if one time series missing, e.g. for:
      %-pilot_RS_cmrr_mbep2d_se_b0_task_12s_off_18s_on_MB1_20200217155721_9
      system("pdfjam diagnostic_plots_time_course_mean_all.pdf diagnostic_plots_time_course_mean_visual.pdf diagnostic_plots_time_course_mean_motor.pdf diagnostic_plots_time_course_mean_somatosensory.pdf --nup 4x1 --landscape --outfile diagnostic_plots_time_courses.pdf");
      system('pdfcrop diagnostic_plots_time_courses.pdf diagnostic_plots_time_courses.pdf');
      
      
      %-analysing volume outliers
      
      load time_course_mean_visual
      
      motion_pars_1 = dir('rp*.txt');
      motion_pars_1 = motion_pars_1(1).name;
      motion_pars_1 = readtable(motion_pars_1);
      motion_pars_1 = table2array(motion_pars_1);
      
      if make_adc == 1
         motion_pars_2 = dir('rp*.txt');
         motion_pars_2 = motion_pars_2(2).name;
         motion_pars_2 = readtable(motion_pars_2);
         motion_pars_2 = table2array(motion_pars_2);
      else
         motion_pars_2 = NaN(length(motion_pars_1), 6);
      end
      
      [B, I] = sort(time_course_mean);
      outliers_indices = [I(1:5) I(end-4:end)];
      
      figure('visible', 'off');
      subplot(3,2,1)
      plot(time_course_mean)
      for i = 1:10; xline(outliers_indices (i), '--', 'LineWidth', 0.01); end
      title('time course mean visual');
      for motion_par_id = 2:6
         subplot(3,2,motion_par_id)
         plot(motion_pars_1(:,motion_par_id)); hold on;
         plot(motion_pars_2(:,motion_par_id))
         for i = 1:10; xline(outliers_indices (i), '--', 'LineWidth', 0.01); end
         title(['motion parameter ' num2str(motion_par_id)]);
      end
      sgtitle(fig_title_text, 'FontSize', 8);
      if analysis_id < 10
         print_to_svg_to_pdf(['outliers_check_0' num2str(analysis_id)]);
      else
         print_to_svg_to_pdf(['outliers_check_'  num2str(analysis_id)]);
      end
      
      %-to save some space
      delete(fullfile([replace_data '.feat'], 'stats', 'res4d.nii.gz'));
      delete(fullfile('SPM_task_GLM_MNI_space',        'res4d.nii.gz'));
      delete(fullfile('SPM_task_GLM_subject_space',    'res4d.nii.gz'));
      
   end
   
   %-to save some space
   
   niftis = dir('**/*.nii');
   for nifti_id = 1:length(niftis)
      system(['gzip ' fullfile(niftis(nifti_id).folder, niftis(nifti_id).name)]);
   end
   
   if exist( 'dfMRI_denoised_unringed_diff.nii.gz',   'file') == 2
      delete('dfMRI_denoised_unringed_diff.nii.gz');
   end
   if exist( 'dfMRI_denoised_unringed_padded.nii.gz', 'file') == 2
      delete('dfMRI_denoised_unringed_padded.nii.gz');
   end
   
   clear func_map
   
end
