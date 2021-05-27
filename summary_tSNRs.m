

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Summarising tSNR maps.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       April 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd /home/wiktor/Desktop/Dfmri/pipeline/analyses_output/


for RS_yn = [1 0]
   
   if RS_yn == 1
      disp(' '); disp('RESTING-STATE    2mm   iso'); disp(' ');
   else
      disp(' '); disp('TASK             2.5mm iso'); disp(' ');
   end
   
   for group_id = 1:4
      
      disp(' ');
      disp(' ');
      
      if group_id == 1
         folders = dir('*_3T*_se*');
         disp('BOLD: 3T');
      elseif group_id == 2
         folders = dir('*_7T*_se*');
         disp('BOLD: 7T');
      elseif group_id == 3
         folders = dir('*_3T*_diff*0pt2_1*');
         disp('ADC: 3T');
      else
         folders = dir('*_7T*_diff*0pt2_1*');
         disp('ADC: 7T');
      end
      
      if RS_yn == 1
         folders  = folders( contains({folders.name}, '_rs_'));
         %-for this subject, RS space confused with task space...
         folders  = folders(~contains({folders.name}, 'BB_004_3T'));
         vox_size = 2;
      else
         folders  = folders( contains({folders.name}, '_18s_off_12s_on_'));
         vox_size = 2.5;
      end
      
      folders = folders(~contains({folders.name}, '_monopolar_'));
      folders = folders(~contains({folders.name}, 'TR350'));
      folders = folders(~contains({folders.name}, 'BB_010_7T'));
      folders = folders(~contains({folders.name}, 'MB1'));
      folders = folders(~contains({folders.name}, 'BB_003_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200710180104_11'));
      folders = folders(~contains({folders.name}, 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200706151255_16'));
      folders = folders(~contains({folders.name}, 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200804162733_15'));
      folders = folders(~contains({folders.name}, 'BB_014_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201029092350_14_combined'));
      folders = folders(~contains({folders.name}, 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201030143609_9_combined'));
      
      ok_vec  = zeros(1, length(folders));
      for folder_id = 1:length(folders)
         if exist(fullfile(folders(folder_id).name, 'ROI_means.mat'), 'file') == 2
            ok_vec(folder_id) = 1;
         end
      end
      folders = folders(ok_vec==1);
      
      disp(' ');
      disp(['No. of subjects: ' num2str(length(folders))]);
      
      tSNR_mean   = NaN(100, 4);
      tSNR_median = NaN(100, 4);
      
      for folder_id = 1:length(folders)
         
         cd(folders(folder_id).name);
         
         brain_mask = niftiread('brain_mask.nii.gz');
         
         header = niftiinfo('SNR/dfMRI_raw_hpf.nii.gz');
         if abs(header.PixelDimensions(1) - vox_size) > 0.1
            cd ..
            continue
         end
         
         if group_id == 1 || group_id == 2
            
            tSNR_before_PCA_denoising = niftiread('SNR/tSNR_before_PCA_denoising.nii.gz');
            tSNR_after_PCA_denoising  = niftiread('SNR/tSNR_after_PCA_denoising.nii.gz');
            
            tSNR_mean  (folder_id, 1) = mean  (tSNR_before_PCA_denoising(brain_mask==1));
            tSNR_mean  (folder_id, 2) = mean  (tSNR_after_PCA_denoising (brain_mask==1));
            
            tSNR_median(folder_id, 1) = median(tSNR_before_PCA_denoising(brain_mask==1));
            tSNR_median(folder_id, 2) = median(tSNR_after_PCA_denoising (brain_mask==1));
            
         else
            
            tSNR_before_PCA_denoising_b_0 = niftiread('SNR/tSNR_before_PCA_denoising_b_0.nii.gz');
            tSNR_after_PCA_denoising_b_0  = niftiread('SNR/tSNR_after_PCA_denoising_b_0.nii.gz');
            tSNR_before_PCA_denoising_b_1 = niftiread('SNR/tSNR_before_PCA_denoising_b_1.nii.gz');
            tSNR_after_PCA_denoising_b_1  = niftiread('SNR/tSNR_after_PCA_denoising_b_1.nii.gz');
            
            tSNR_mean  (folder_id, 1) = mean  (tSNR_before_PCA_denoising_b_0(brain_mask==1));
            tSNR_mean  (folder_id, 2) = mean  (tSNR_after_PCA_denoising_b_0 (brain_mask==1));
            tSNR_mean  (folder_id, 3) = mean  (tSNR_before_PCA_denoising_b_1(brain_mask==1));
            tSNR_mean  (folder_id, 4) = mean  (tSNR_after_PCA_denoising_b_1 (brain_mask==1));
            
            tSNR_median(folder_id, 1) = median(tSNR_before_PCA_denoising_b_0(brain_mask==1));
            tSNR_median(folder_id, 2) = median(tSNR_after_PCA_denoising_b_0 (brain_mask==1));
            tSNR_median(folder_id, 3) = median(tSNR_before_PCA_denoising_b_1(brain_mask==1));
            tSNR_median(folder_id, 4) = median(tSNR_after_PCA_denoising_b_1 (brain_mask==1));
            
         end
         
         cd ..
         
      end
      
      disp(['Lower  b-value before MP-PCA denoising mean: ' num2str(round(mean(tSNR_mean(:,1), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_mean(:,1), 'omitnan'))) ',    median: ' num2str(round(mean(tSNR_median(:,1), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_median(:,1), 'omitnan')))]);
      disp(['Lower  b-value after  MP-PCA denoising mean: ' num2str(round(mean(tSNR_mean(:,2), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_mean(:,2), 'omitnan'))) ',    median: ' num2str(round(mean(tSNR_median(:,2), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_median(:,2), 'omitnan')))]);
      
      if group_id > 2
         disp(['Higher b-value before MP-PCA denoising mean: ' num2str(round(mean(tSNR_mean(:,3), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_mean(:,3), 'omitnan'))) ',    median: ' num2str(round(mean(tSNR_median(:,3), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_median(:,3), 'omitnan')))]);
         disp(['Higher b-value after  MP-PCA denoising mean: ' num2str(round(mean(tSNR_mean(:,4), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_mean(:,4), 'omitnan'))) ',    median: ' num2str(round(mean(tSNR_median(:,4), 'omitnan'))) '  sd: ' num2str(round(std(tSNR_median(:,4), 'omitnan')))]);
      end
      
   end
   
end
