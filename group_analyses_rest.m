

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Performing resting-state (RS) group analyses: plotting group-averaged functional
%%%%   connectivity (FC) matrices, as well as running statistical tests. Run this code after
%%%%   performing single-subject analyses with 'fMRI_processing_single_subject.m'.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       October 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage_software = '/home/wiktor/Desktop/Dfmri/pipeline/software';
path_output          = '/home/wiktor/Desktop/Dfmri/pipeline/analyses_output';


for clean_yn = 0:1
   
   if clean_yn == 0
      clean_suffix = '';
   else
      clean_suffix = '_clean';
   end
   
   cd(path_output);
   
   folders = dir('BB*');
   
   imgs_se_3T               = string([]);
   imgs_se_7T               = string([]);
   imgs_diff_0_1_3T         = string([]);
   imgs_diff_0_1_7T         = string([]);
   imgs_diff_0pt29_1pt45_3T = string([]);
   imgs_diff_0pt29_1pt45_7T = string([]);
   imgs_diff_0pt2_1_mono_3T = string([]);
   imgs_diff_0pt2_1_mono_7T = string([]);
   imgs_diff_0pt2_1_bi_3T   = string([]);
   imgs_diff_0pt2_1_bi_7T   = string([]);
   
   for folder_id = 1:length(folders)
      
      img = folders(folder_id).name;
      
      %-this analysis is only for rest data
      if contains(img, '_task_')
         continue
      end
      
      %-for this subject, RS space confused with task space...
      if contains(img, 'BB_004_3T')
         continue
      end
      
      %-identified by Yujian as bad image
      if contains(img, 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_1_rs_bipolar_20210219110948_13')
         continue
      end
      
      if isempty(dir(fullfile(img, ['ROI_means' clean_suffix '_GSR.mat'])))
         disp(['problem with ' img]);
         continue
      end
      
      if contains(img, '_se_')
         
         if     contains(img, '_3T_')
            imgs_se_3T(end+1) = img;
         elseif contains(img, '_7T_')
            imgs_se_7T(end+1) = img;
         end
         
      elseif contains(img, '_diff_')
         
         %-no RS b0_1 at 3T
         if     contains(img, 'BB_001_7T') || contains(img, 'BB_002_7T')
            imgs_diff_0_1_7T(end+1)         = img;
         elseif contains(img, '_3T_')      && contains(img, '0pt29_1pt45')
            imgs_diff_0pt29_1pt45_3T(end+1) = img;
         elseif contains(img, '_7T_')      && contains(img, '0pt29_1pt45')
            imgs_diff_0pt29_1pt45_7T(end+1) = img;
         elseif contains(img, '_3T_')      && contains(img, '_monopolar_')
            imgs_diff_0pt2_1_mono_3T(end+1) = img;
         elseif contains(img, '_7T_')      && contains(img, '_monopolar_')
            imgs_diff_0pt2_1_mono_7T(end+1) = img;
         elseif contains(img, '_3T_')
            imgs_diff_0pt2_1_bi_3T(end+1)   = img;
         elseif contains(img, '_7T_')
            imgs_diff_0pt2_1_bi_7T(end+1)   = img;
         end
         
      end
      
   end
   
   imgs_se_3T               = imgs_se_3T';
   imgs_se_7T               = imgs_se_7T';
   imgs_diff_0_1_3T         = imgs_diff_0_1_3T';
   imgs_diff_0_1_7T         = imgs_diff_0_1_7T';
   imgs_diff_0pt29_1pt45_3T = imgs_diff_0pt29_1pt45_3T';
   imgs_diff_0pt29_1pt45_7T = imgs_diff_0pt29_1pt45_7T';
   imgs_diff_0pt2_1_mono_3T = imgs_diff_0pt2_1_mono_3T';
   imgs_diff_0pt2_1_mono_7T = imgs_diff_0pt2_1_mono_7T';
   imgs_diff_0pt2_1_bi_3T   = imgs_diff_0pt2_1_bi_3T';
   imgs_diff_0pt2_1_bi_7T   = imgs_diff_0pt2_1_bi_7T';
   
   fileID = fopen(fullfile(path_manage_software, 'brain_parcellation', 'labels_Neuromorphometrics.txt'));
   A      = textscan(fileID, '%s', 'Delimiter', '', 'headerLines', 1);
   A      = A{1};
   
   ROIs = zeros (1, size(A, 1));
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
   
   
   ROIs_intersection = 1:length(ROIs);
   
   for type_id = 1:8
      
      if     type_id == 1
         imgs = imgs_se_3T;
      elseif type_id == 2
         imgs = imgs_se_7T;
      elseif type_id == 3
         imgs = imgs_diff_0_1_3T;
      elseif type_id == 4
         imgs = imgs_diff_0_1_7T;
      elseif type_id == 5
         imgs = imgs_diff_0pt2_1_bi_3T;
      elseif type_id == 6
         imgs = imgs_diff_0pt2_1_bi_7T;
      elseif type_id == 7
         imgs = imgs_diff_0pt2_1_mono_3T;
      elseif type_id == 8
         imgs = imgs_diff_0pt2_1_mono_7T;
      end
      for img_id = 1:length(imgs)
         cd(imgs(img_id));
         load ROI_means_GSR
         cd ..
         ROIs_intersection = intersect(ROIs_intersection, find(ROI_means(1,:)~=0 & ~isnan(ROI_means(1,:))));
      end
      
   end
   
   
   for cleaning_with_GSR_yn = 0:1
      
      if clean_yn == 0 && cleaning_with_GSR_yn == 0
         continue
      end
      
      for map_0pt2_yn = 0:1
         
         if map_0pt2_yn == 0
            map_0pt2_suffix = '';
         else
            map_0pt2_suffix = '_map_0pt2';
         end
         
         for bi_mono_id = 1:2
            
            figure('visible', 'off', 'pos', [0 0 600 400]);
            
            ROI_group_corrs_non_empty_merged = {};
            
            for subplot_id = 1:4
               
               if subplot_id == 1
                  subplot(2, 3, 1);
                  if map_0pt2_yn == 0
                     imgs = imgs_se_3T;
                  else
                     if bi_mono_id == 1
                        imgs = imgs_diff_0pt2_1_bi_3T;
                     else
                        imgs = imgs_diff_0pt2_1_mono_3T;
                     end
                  end
                  if map_0pt2_yn == 0
                     title_t = 'SE BOLD 3T';
                  else
                     title_t = 'pseudo SE BOLD 3T';
                  end
               elseif subplot_id == 2
                  subplot(2, 3, 2);
                  if bi_mono_id == 1
                     imgs = imgs_diff_0pt2_1_bi_3T;
                  else
                     imgs = imgs_diff_0pt2_1_mono_3T;
                  end
                  title_t = 'ADC 0.2/1 3T';
               elseif subplot_id == 3
                  subplot(2, 3, 4);
                  if map_0pt2_yn == 0
                     imgs = imgs_se_7T;
                  else
                     if bi_mono_id == 1
                        imgs = imgs_diff_0pt2_1_bi_7T;
                     else
                        imgs = imgs_diff_0pt2_1_mono_7T;
                     end
                  end
                  if map_0pt2_yn == 0
                     title_t = 'SE BOLD 7T';
                  else
                     title_t = 'pseudo SE BOLD 7T';
                  end
               elseif subplot_id == 4
                  subplot(2, 3, 5);
                  if bi_mono_id == 1
                     imgs = imgs_diff_0pt2_1_bi_7T;
                  else
                     imgs = imgs_diff_0pt2_1_mono_7T;
                  end
                  title_t = 'ADC 0.2/1 7T';
               end
               
               if     bi_mono_id == 1 && (subplot_id == 2 || subplot_id == 4)
                  bi_mono_suffix = ' bipolar';
               elseif bi_mono_id == 2 && (subplot_id == 2 || subplot_id == 4)
                  bi_mono_suffix = ' monopolar';
               else
                  bi_mono_suffix = '';
               end
               
               ROI_group_corrs      = zeros(length(ROIs_t), length(ROIs_t));
               No_subjects_per_corr = zeros(length(ROIs_t), length(ROIs_t));
               No_of_missing        = 0;
               
               for img_id = 1:length(imgs)
                  
                  cd(imgs(img_id));
                  %-here we identify the imgs with too many volume outliers, so that the imgs were not processed
                  if exist(['ROI_means' clean_suffix '_GSR.mat'], 'file') ~= 2
                     disp(['no ROI_means.mat for ' imgs(img_id)]);
                     No_of_missing = No_of_missing + 1;
                     cd ..
                     continue;
                  end
                  
                  if subplot_id == 1 || subplot_id == 3
                     if cleaning_with_GSR_yn == 1
                        load(['ROI_means' clean_suffix map_0pt2_suffix '_GSR']);
                     else
                        load(['ROI_means' clean_suffix map_0pt2_suffix]);
                     end
                  else
                     if cleaning_with_GSR_yn == 1
                        load(['ROI_means' clean_suffix '_GSR']);
                     else
                        load(['ROI_means' clean_suffix]);
                     end
                  end
                  
                  cd ..
                  ROI_means(isnan(ROI_means)) = 0;
                  corr_ROI_means = corr(ROI_means);
                  %-to include ROI pairs where less than all subjects available, uncomment:
                  %-corr_ROI_means(isnan(corr_ROI_means)) = 0;
                  ROI_group_corrs = ROI_group_corrs + corr_ROI_means;
                  No_subjects_per_corr(corr(ROI_means)~=0 & ~isnan(corr(ROI_means))) = No_subjects_per_corr(corr(ROI_means)~=0 & ~isnan(corr(ROI_means))) + 1;
                  
               end
               
               non_empty_indices = ROIs_intersection;
               ROI_group_corrs_non_empty = ROI_group_corrs(non_empty_indices, non_empty_indices) ./ No_subjects_per_corr(non_empty_indices, non_empty_indices);
               imagesc(ROI_group_corrs_non_empty);
               vec_aux                                      = reshape(ROI_group_corrs_non_empty, 1, []);
               ROI_group_corrs_non_empty_merged{subplot_id} = vec_aux(reshape(tril(ones(length(ROIs_intersection)), -1), 1, [])==1);
               caxis([-1 1]);
               colormap('jet');
               
               ax            = gca;
               ax.YTick      = 1:length(non_empty_indices);
               ax.XTick      = 1:length(non_empty_indices);
               %-the spaces in YTickLabels appear when exporting to pdf, it is a known problem in matlab
               %-I can get rid of it changing the interpreter, but then Helvetica font doesn't work
               % ax.YTickLabel = ROIs_t(non_empty_indices);
               % ax.XTickLabel = ROIs_t(non_empty_indices);
               ax.YTickLabel = repmat({'                  '}, 1, length(non_empty_indices));
               ax.XTickLabel = repmat({'                  '}, 1, length(non_empty_indices));
               xtickangle(90);
               set(ax, 'xaxisLocation', 'top');
               set(ax, 'FontSize', 1.8);
               
               if contains(title_t, '7T')
                  plot_pos    = get(gca, 'position');
                  plot_pos(2) = plot_pos(2)+0.02;
                  set(gca, 'position', plot_pos);
               end
               
               if subplot_id == 2 || subplot_id == 4
                  t = title({title_t(end-2:end); ' '; [title_t(1:end-3) bi_mono_suffix ' (' num2str(length(imgs)-No_of_missing) ' subjects)']}, 'FontSize', 6, 'interpreter', 'latex');
               else
                  t = title({' ';                ' '; [title_t(1:end-3) bi_mono_suffix ' (' num2str(length(imgs)-No_of_missing) ' subjects)']}, 'FontSize', 6, 'interpreter', 'latex');
               end
               title_pos    = get(t, 'position');
               title_pos(2) = title_pos(2)+3;
               set(t, 'position', title_pos);
               
               if subplot_id == 3
                  cb = colorbar;
                  cb.Position(1) = 0.36; % or + [left, bottom, width, height] to place it where you want
                  set(cb, 'FontSize', 5);
               end
               
               if subplot_id == 2 || subplot_id == 4
                  
                  subplot(2,3,1.5*subplot_id);
                  conn_SE_BOLD = ROI_group_corrs_non_empty_merged{subplot_id-1};
                  conn_dfMRI   = ROI_group_corrs_non_empty_merged{subplot_id};
                  scatter(conn_SE_BOLD, conn_dfMRI, 2, 'filled');
                  set(gca, 'FontSize', 4);
                  if contains(title_t, '7T')
                     plot_pos    = get(gca, 'position');
                     plot_pos(2) = plot_pos(2)+0.02;
                     set(gca, 'position', plot_pos);
                  end
                  
                  if map_0pt2_yn == 0
                     xlabel('SE BOLD connectivity',        'interpreter', 'latex');
                  else
                     xlabel('pseudo SE BOLD connectivity', 'interpreter', 'latex');
                  end
                  
                  ylabel('ADC 0.2/1 connectivity', 'interpreter', 'latex');
                  xlim([-1.05 1.05]);
                  ylim([-1.05 1.05]);
                  xline(0);
                  yline(0);
                  
                  which_negative = find(conn_SE_BOLD < 0);
                  which_positive = find(conn_SE_BOLD > 0);
                  LM_neg = fitlm(conn_SE_BOLD(which_negative), conn_dfMRI(which_negative));
                  LM_pos = fitlm(conn_SE_BOLD(which_positive), conn_dfMRI(which_positive));
                  hold on; plot(conn_SE_BOLD(which_negative), predict(LM_neg, conn_SE_BOLD(which_negative)'), 'red', 'LineWidth', 1.25);
                  hold on; plot(conn_SE_BOLD(which_positive), predict(LM_pos, conn_SE_BOLD(which_positive)'), 'red', 'LineWidth', 1.25);
                  text(-0.7, 0.86, ['slope=' num2str(round(LM_neg.Coefficients.Estimate(2), 2))], 'FontSize', 6, 'interpreter', 'latex');
                  text( 0.1, 0.86, ['slope=' num2str(round(LM_pos.Coefficients.Estimate(2), 2))], 'FontSize', 6, 'interpreter', 'latex');
                  set(gca, 'FontSize', 5);
                  box on;
                  
               end
               
               pause(5);
               
            end
            
            cd ../figures
            if bi_mono_id == 1
               print([         'group_rest_correlations_abstract_bipolar_GSR_'   num2str(cleaning_with_GSR_yn) clean_suffix map_0pt2_suffix], '-dpdf');
               pause(5);
               system(['pdfcrop group_rest_correlations_abstract_bipolar_GSR_'   num2str(cleaning_with_GSR_yn) clean_suffix map_0pt2_suffix '.pdf group_rest_correlations_abstract_bipolar_GSR_'   num2str(cleaning_with_GSR_yn) clean_suffix map_0pt2_suffix '.pdf']);
            else
               print([         'group_rest_correlations_abstract_monopolar_GSR_' num2str(cleaning_with_GSR_yn) clean_suffix map_0pt2_suffix], '-dpdf');
               pause(5);
               system(['pdfcrop group_rest_correlations_abstract_monopolar_GSR_' num2str(cleaning_with_GSR_yn) clean_suffix map_0pt2_suffix '.pdf group_rest_correlations_abstract_monopolar_GSR_' num2str(cleaning_with_GSR_yn) clean_suffix map_0pt2_suffix '.pdf']);
            end
            
            pause(5);
            cd ../analyses_output
            
         end
         
      end
      
   end
   
   
   if clean_yn == 0
      
      
      for field_strength = [3 7]
         
         figure('visible', 'off', 'pos', [0 0 600 600]);
         
         for subplot_id = 1:9
            
            subplot(3, 3, subplot_id);
            
            title_t = '';
            y_label = '';
            
            if subplot_id == 1
               y_label = 'SE BOLD';
               title_t = 'raw';
               imgs    = [imgs_se_3T; imgs_se_7T];
               suffix  = '';
            elseif subplot_id == 2
               title_t = 'GSR';
               imgs    = [imgs_se_3T; imgs_se_7T];
               suffix  = '_GSR';
            elseif subplot_id == 3
               title_t = 'ICA cleaning + GSR';
               imgs    = [imgs_se_3T; imgs_se_7T];
               suffix  = '_clean_GSR';
            elseif subplot_id == 4
               y_label = 'pseudo SE BOLD';
               imgs    = [imgs_diff_0pt2_1_bi_3T; imgs_diff_0pt2_1_bi_7T];
               suffix  = '_map_0pt2';
            elseif subplot_id == 5
               imgs    = [imgs_diff_0pt2_1_bi_3T; imgs_diff_0pt2_1_bi_7T];
               suffix  = '_map_0pt2_GSR';
            elseif subplot_id == 6
               imgs    = [imgs_diff_0pt2_1_bi_3T; imgs_diff_0pt2_1_bi_7T];
               suffix  = '_clean_map_0pt2_GSR';
            elseif subplot_id == 7
               y_label = 'ADC 0.2/1 bipolar';
               imgs    = [imgs_diff_0pt2_1_bi_3T; imgs_diff_0pt2_1_bi_7T];
               suffix  = '';
            elseif subplot_id == 8
               imgs    = [imgs_diff_0pt2_1_bi_3T; imgs_diff_0pt2_1_bi_7T];
               suffix  = '_GSR';
            elseif subplot_id == 9
               imgs    = [imgs_diff_0pt2_1_bi_3T; imgs_diff_0pt2_1_bi_7T];
               suffix  = '_clean_GSR';
            end
            
            imgs = imgs(contains(imgs, ['_' num2str(field_strength) 'T_']));
            
            ROI_group_corrs      = zeros(length(ROIs_t), length(ROIs_t));
            No_subjects_per_corr = zeros(length(ROIs_t), length(ROIs_t));
            No_of_missing        = 0;
            
            for img_id = 1:length(imgs)
               
               cd(imgs(img_id));
               %-here we identify the imgs with too many volume outliers, so that the imgs were not processed
               if exist(['ROI_means' suffix '.mat'], 'file') ~= 2
                  disp(['no ROI_means' suffix '.mat for ' char(imgs(img_id))]);
                  No_of_missing = No_of_missing + 1;
                  cd ..
                  continue;
               end
               
               load(['ROI_means' suffix]);
               cd ..
               ROI_means(isnan(ROI_means)) = 0;
               corr_ROI_means = corr(ROI_means);
               ROI_group_corrs = ROI_group_corrs + corr_ROI_means;
               No_subjects_per_corr(corr(ROI_means)~=0 & ~isnan(corr(ROI_means))) = No_subjects_per_corr(corr(ROI_means)~=0 & ~isnan(corr(ROI_means))) + 1;
               
            end
            
            non_empty_indices = ROIs_intersection;
            ROI_group_corrs_non_empty = ROI_group_corrs(non_empty_indices, non_empty_indices) ./ No_subjects_per_corr(non_empty_indices, non_empty_indices);
            imagesc(ROI_group_corrs_non_empty);
            caxis([-1 1]);
            colormap('jet');
            
            ax            = gca;
            ax.YTick      = 1:length(non_empty_indices);
            ax.XTick      = 1:length(non_empty_indices);
            ax.YTickLabel = repmat({'                  '}, 1, length(non_empty_indices));
            ax.XTickLabel = repmat({'                  '}, 1, length(non_empty_indices));
            xtickangle(90);
            set(ax, 'xaxisLocation', 'top');
            set(ax, 'FontSize', 1.8);
            plot_pos    = get(gca, 'position');
            plot_pos(2) = plot_pos(2) + floor((subplot_id-1)/3)*0.03;
            set(gca, 'position', plot_pos);
            
            if subplot_id == 2 || subplot_id == 5 || subplot_id == 8
               No_t = [num2str(length(imgs)-No_of_missing) ' subjects'];
            else
               No_t = '';
            end
            
            t = title({title_t; ''; No_t}, 'FontSize', 8, 'interpreter', 'latex');
            title_pos    = get(t, 'position');
            title_pos(2) = title_pos(2)+3;
            set(t, 'position', title_pos);
            
            ylabel(y_label, 'FontSize', 8, 'interpreter', 'latex');
            
            if subplot_id == 4
               cb = colorbar;
               cb.Position(1) = 0.36;
               set(cb, 'FontSize', 6);
            end
            
         end
         
         cd ../figures
         
         print(         ['group_rest_correlations_3x3_' num2str(field_strength) 'T'], '-dpdf');
         pause(5);
         system(['pdfcrop group_rest_correlations_3x3_' num2str(field_strength) 'T.pdf group_rest_correlations_3x3_' num2str(field_strength) 'T.pdf']);
         pause(5);
         cd ../analyses_output
         
      end
      
   end
   
   
   if false
      
      %-2-way ANOVA
      %-results difficult to interpret -> t-tests more straightforward
      
      No_reps           = 10;
      all_anova2_p_vals = NaN(length(non_empty_indices), length(non_empty_indices), 3);
      cd(path_output);
      
      for ROI_1_id = 1:(length(non_empty_indices)-1)
         disp(ROI_1_id);
         for ROI_2_id = (ROI_1_id + 1) : length(non_empty_indices)
            matrix_for_anova2 = NaN(No_reps*2, 2);
            for row_id = 1:2
               for column_id = 1:2
                  if     row_id == 1 && column_id == 1
                     imgs = imgs_se_3T;
                  elseif row_id == 1 && column_id == 2
                     imgs = imgs_diff_0pt2_1_bi_3T;
                  elseif row_id == 2 && column_id == 1
                     imgs = imgs_se_7T;
                  elseif row_id == 2 && column_id == 2
                     imgs = imgs_diff_0pt2_1_bi_7T;
                  end
                  %-for balanced 2-way ANOVA we need samples of the same size
                  for rep_id = 1:No_reps
                     cd(imgs(rep_id));
                     if (row_id == 1 && column_id == 1) || (row_id == 2 && column_id == 1)
                        load(['ROI_means' clean_suffix map_0pt2_suffix '_GSR']);
                     else
                        load(['ROI_means' clean_suffix '_GSR']);
                     end
                     cd ..
                     ROI_means(isnan(ROI_means)) = 0;
                     corr_ROI_means = corr(ROI_means);
                     matrix_for_anova2((row_id-1)*No_reps + rep_id, column_id) = corr_ROI_means(non_empty_indices(ROI_1_id), non_empty_indices(ROI_2_id));
                  end
               end
            end
            all_anova2_p_vals(ROI_1_id, ROI_2_id, :) = anova2(matrix_for_anova2, No_reps, 'off');
         end
      end
      
      No_obs_with_low_p = 40;
      ROIs_t_non_empty  = ROIs_t(non_empty_indices);
      
      for level_id = 1:3
         
         if level_id == 1
            fprintf('\n\n\nfield difference\n\n\n');
         elseif level_id == 2
            fprintf('\n\n\nmodality difference\n\n\n');
         elseif level_id == 3
            fprintf('\n\n\ninteraction\n\n\n');
         end
         [B, I] = sort(reshape(all_anova2_p_vals(:,:,level_id), 1, []), 'ascend');
         for obs_id = 1:No_obs_with_low_p
            ROI_1_id = fix((I(obs_id)-1)/length(non_empty_indices)) + 1;
            ROI_2_id = I(obs_id) - (ROI_1_id-1)*length(non_empty_indices);
            %-with Bonferroni correction
            disp([char(sprintf("%0.5f", 1035*B(obs_id))) '    ' char(sprintf("%-43s", char(ROIs_t_non_empty(ROI_1_id)))) ' vs. ' char(ROIs_t_non_empty(ROI_2_id))]);
         end
         
      end
      
   end
   
   
   %-t-tests
   
   non_empty_indices = ROIs_intersection;
   all_t_test_p_vals = NaN(length(non_empty_indices), length(non_empty_indices), 4);
   No_obs_with_low_p = 20;
   ROIs_t_non_empty  = ROIs_t(non_empty_indices);
   
   for case_id = 1:4
      
      if     case_id == 1
         imgs_1 = imgs_se_3T;
         imgs_2 = imgs_se_7T;
         case_t = 'SE BOLD: 3T/7T';
      elseif case_id == 2
         imgs_1 = imgs_diff_0pt2_1_bi_3T;
         imgs_2 = imgs_diff_0pt2_1_bi_7T;
         case_t = 'ADC: 3T/7T';
      elseif case_id == 3
         imgs_1 = imgs_se_3T;
         imgs_2 = imgs_diff_0pt2_1_bi_3T;
         case_t = '3T: SE BOLD/ADC';
      elseif case_id == 4
         imgs_1 = imgs_se_7T;
         imgs_2 = imgs_diff_0pt2_1_bi_7T;
         case_t = '7T: SE BOLD/ADC';
      end
      
      all_corrs_1  = NaN(length(non_empty_indices), length(non_empty_indices), length(imgs_1));
      all_corrs_2  = NaN(length(non_empty_indices), length(non_empty_indices), length(imgs_2));
      
      for img_1_id = 1:length(imgs_1)
         cd(imgs_1(img_1_id));
         load(['ROI_means' clean_suffix '_GSR']);
         cd ..
         ROI_means = ROI_means(:,non_empty_indices);
         %-with Fisher transformation
         %-https://blogs.sas.com/content/iml/2017/09/20/fishers-transformation-correlation.html
         all_corrs_1(:, :, img_1_id) = atanh(corr(ROI_means));
      end
      
      for img_2_id = 1:length(imgs_2)
         cd(imgs_2(img_2_id));
         load(['ROI_means' clean_suffix '_GSR']);
         cd ..
         ROI_means = ROI_means(:,non_empty_indices);
         %-with Fisher transformation
         all_corrs_2(:, :, img_2_id) = atanh(corr(ROI_means));
      end
      
      for ROI_1_id = 1:(length(non_empty_indices)-1)
         for ROI_2_id = (ROI_1_id+1):length(non_empty_indices)
            [h, p] = ttest2(all_corrs_1(ROI_1_id, ROI_2_id, :), all_corrs_2(ROI_1_id, ROI_2_id, :));
            all_t_test_p_vals(ROI_1_id, ROI_2_id, case_id) = p;
         end
      end
      
      disp(' ');
      disp(case_t);
      
      [B, I] = sort(reshape(all_t_test_p_vals(:,:,case_id), 1, []), 'ascend');
      
      for obs_id = 1:No_obs_with_low_p
         
         ROI_1_id = fix((I(obs_id)-1)/length(non_empty_indices)) + 1;
         ROI_2_id = I(obs_id) - (ROI_1_id-1)*length(non_empty_indices);
         %-with Bonferroni correction, 46*(46-1)/2 = 1035; there are 46 ROIs
         %-the first way of outputting (commented out) is good for looking at the results in the terminal
         %-the second way of outputting (not commented out) is good for exporting to Excel ('!' can separate columns)
         %{
			out_1 = char(sprintf("%-43s", char(ROIs_t_non_empty(ROI_1_id))));
			out_2 = char(sprintf("%-43s", char(ROIs_t_non_empty(ROI_2_id))));
			out_3 = char(sprintf("%0.3f",  mean(all_corrs_1(ROI_1_id, ROI_2_id, :))));
			out_4 = char(sprintf("%0.3f",  mean(all_corrs_2(ROI_1_id, ROI_2_id, :))));
         %}
         out_1 = ['!' char(ROIs_t_non_empty(ROI_1_id))];
         out_2 = ['!' char(ROIs_t_non_empty(ROI_2_id))];
         out_3 = ['!' char(sprintf("%0.3f",  mean(all_corrs_1(ROI_1_id, ROI_2_id, :))))];
         out_4 = ['!' char(sprintf("%0.3f",  mean(all_corrs_2(ROI_1_id, ROI_2_id, :))))];
         if out_3(1) ~= '-'
            out_3 = [' ' out_3];
         end
         if out_4(1) ~= '-'
            out_4 = [' ' out_4];
         end
         disp([char(sprintf("%0.5f", 1035*B(obs_id))) '   ' out_1  ' vs. ' out_2 '  ' out_3 '  ' out_4 ]);
         
      end
      
   end
   
end

cd ..
