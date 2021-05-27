

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Performing spatial task group analyses: averaging z-statistic maps across subjects.
%%%%   Looking how much activity there is in different ROIs. These ROIs are from the NMM atlas
%%%%   and the JHU WM atlas. Run this code after performing single-subject analyses with
%%%%   'fMRI_processing_single_subject.m', and after running simplified GLM analyses with
%%%%   'GLM_single_subject_one_basis_function_only.m'.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       January 2021 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage          = '/home/wiktor/Desktop/Dfmri/pipeline';
path_manage_software = fullfile(path_manage, 'software');


%-loading SPM
addpath(genpath(fullfile(path_manage_software, 'spm12')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-reading the Neuromorphometrics atlas parcellation

fileID = fopen(fullfile(path_manage_software, 'brain_parcellation', 'labels_Neuromorphometrics_JHU.txt'));
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
ROIs_t = cellstr(ROIs_t);


for analysis_id = 1:3
   
   cd(fullfile(path_manage, 'analyses_output'));
   
   if analysis_id == 1
      
      zfstat_files_all = cellstr(spm_select('FPListRec', fullfile(path_manage, 'analyses_output'), '^zfstat1_MNI.nii.gz'));
      analysis_t       = 'FIR';
      
   elseif analysis_id == 2
      
      zfstat_files_all = cellstr(spm_select('FPListRec', fullfile(path_manage, 'analyses_output'), '^zstat_MNI.nii.gz'));
      zfstat_files_all = zfstat_files_all(contains(zfstat_files_all, 'SPM_task_GLM_subject_space_one_rectangular_basis_function'));
      analysis_t       = 'rectangular_basis_function';
      
   elseif analysis_id == 3
      
      zfstat_files_all = cellstr(spm_select('FPListRec', fullfile(path_manage, 'analyses_output'), '^zstat_MNI.nii.gz'));
      zfstat_files_all = zfstat_files_all(contains(zfstat_files_all, 'SPM_task_GLM_subject_space_one_canonical_basis_function'));
      analysis_t       = 'canonical_basis_function';
      
   end
   
   disp(' '); disp(analysis_t); disp(' ');
   
   zfstat_files_all = zfstat_files_all(contains(zfstat_files_all, '_18s_off_12s_on_'));
   
   %-subject 'BB_010_7T' was partially sleeping, at least at the end (for SE)
   zfstat_files_all = zfstat_files_all(~contains(zfstat_files_all, '_monopolar_') & ~contains(zfstat_files_all, 'TR350') & ~contains(zfstat_files_all, 'BB_010_7T') & ~contains(zfstat_files_all, 'MB1'));
   
   %-additional exclusions (artifacts)
   zfstat_files_all = zfstat_files_all(~contains(zfstat_files_all, 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200706151255_16') & ...
      ~contains(zfstat_files_all, 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200804162733_15') & ...
      ~contains(zfstat_files_all, 'BB_014_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201029092350_14_combined') & ...
      ~contains(zfstat_files_all, 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201030143609_9_combined'));
   
   for case_id = 1:4
      
      if case_id == 1
         case_name = 'SE_BOLD_3T';
      elseif case_id == 2
         case_name = 'ADC_3T';
      elseif case_id == 3
         case_name = 'SE_BOLD_7T';
      elseif case_id == 4
         case_name = 'ADC_7T';
      end
      
      if case_id == 1 || case_id == 2
         zfstat_files = zfstat_files_all(contains(zfstat_files_all, '_3T_'));
      else
         zfstat_files = zfstat_files_all(contains(zfstat_files_all, '_7T_'));
      end
      
      if case_id == 1 || case_id == 3
         zfstat_files = zfstat_files(contains(zfstat_files, '_se_'));
      else
         zfstat_files = zfstat_files(contains(zfstat_files, '0pt2_1'));
      end
      
      if     analysis_id == 1 && (case_id == 1 || case_id == 3)
         zfstat_files = zfstat_files(contains(zfstat_files, 'ardfMRI_denoised_unringed_sdc.feat'));
      elseif analysis_id == 1 && (case_id == 2 || case_id == 4)
         zfstat_files = zfstat_files(contains(zfstat_files, 'adc.feat'));
      elseif analysis_id == 2
         zfstat_files = zfstat_files(contains(zfstat_files, 'SPM_task_GLM_subject_space_one_rectangular_basis_function'));
      elseif analysis_id == 3
         zfstat_files = zfstat_files(contains(zfstat_files, 'SPM_task_GLM_subject_space_one_canonical_basis_function'));
      end
      
      fprintf([case_name '\n']);
      disp(['No. of subjects: ' num2str(length(zfstat_files))]);
      
      No_of_vox_per_ROI     = NaN([length(ROIs) length(zfstat_files)]);
      No_of_act_vox_per_ROI = NaN([length(ROIs) length(zfstat_files)]);
      
      for file_id = 1:length(zfstat_files)
         
         cd(fileparts(zfstat_files{file_id}));
         
         if analysis_id == 1
            zfstat = reshape(niftiread('zfstat1.nii.gz'), 1, []);
            header = niftiinfo('zfstat1.nii.gz');
            cd ../..
         else
            zfstat = reshape(niftiread('zstat.nii.gz'), 1, []);
            header = niftiinfo('zstat.nii.gz');
            cd ..
         end
         
         NMM_atlas  = reshape(niftiread('NMM_JHU_atlas_subject_space.nii.gz'), 1, []);
         for ROI_id = 1:length(ROIs)
            if sum(NMM_atlas == ROIs(ROI_id)) == 0
               continue
            else
               %-some datasets 2mm iso, some 2.5mm iso
               %-we scale numbers of voxels, pretending to have 1mm iso voxels
               No_of_vox_per_ROI    (ROI_id, file_id) = sum(NMM_atlas == ROIs(ROI_id))               * prod(header.PixelDimensions);
               No_of_act_vox_per_ROI(ROI_id, file_id) = sum(zfstat(NMM_atlas == ROIs(ROI_id)) > 3.1) * prod(header.PixelDimensions);
            end
         end
         
      end
      
      %-means calculated only for ROIs which are covered by all datasets
      mean_No_of_vox_per_ROI       = mean(No_of_vox_per_ROI,       2);
      mean_No_of_act_vox_per_ROI   = mean(No_of_act_vox_per_ROI,   2);
      median_No_of_act_vox_per_ROI = median(No_of_act_vox_per_ROI, 2);
      
      [B, I]                       = sort(mean_No_of_act_vox_per_ROI ./ mean_No_of_vox_per_ROI, 'descend', 'MissingPlacement', 'last');
      
      disp('perc     mean     median   ROI');
      disp(' ');
      for ind = 1:length(ROIs)
         
         if isnan(B(ind))
            continue
         end
         perc_aux   = num2str(round(100*B(ind)));
         mean_aux   = num2str(round(mean_No_of_act_vox_per_ROI(I(ind))));
         median_aux = num2str(round(median_No_of_act_vox_per_ROI(I(ind))));
         %-the first way of outputting (commented out) is good for looking at the results in the terminal
         %-the second way of outputting (not commented out) is good for exporting to Excel ('!' can separate columns)
         %{
			row_str    = [        perc_aux   repmat(' ', 1, 9-length(perc_aux))];
			row_str    = [row_str mean_aux   repmat(' ', 1, 9-length(mean_aux))];
			row_str    = [row_str median_aux repmat(' ', 1, 9-length(median_aux))];
			row_str    = [row_str char(ROIs_t(I(ind)))];
         %}
         row_str    = [        perc_aux   '!'];
         row_str    = [row_str mean_aux   '!'];
         row_str    = [row_str median_aux '!'];
         row_str    = [row_str char(ROIs_t(I(ind)))];
         disp(row_str);
         
      end
      
      fprintf('\n\n\n');
      
      cd(fullfile(path_manage, 'group_whole_brain_analyses'));
      mkdir(analysis_t);
      cd(analysis_t);
      
      system(['fslmerge -t zfstat_' case_name '_merged ' strjoin(zfstat_files)]);
      system(['fslmaths zfstat_' case_name '_merged -Tmean zfstat_' case_name '_mean']);
      system(['fslmaths zfstat_' case_name '_mean -mas /usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz zfstat_' case_name '_mean_masked']);
      
   end
   
end
