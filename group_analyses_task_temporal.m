

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Performing temporal task group analyses: plotting response functions. This analysis is
%%%%   run for voxels selected based on their z-statistics (transformed from F-statistics).
%%%%   By default, cluster thresholding is not performed. To apply the cluster threshold, change
%%%%   'cluster_suffix' to '_clustered'. Run this code after performing single-subject analyses
%%%%   with 'fMRI_processing_single_subject.m'.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       August 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage    = '/home/wiktor/Desktop/Dfmri/pipeline';
folders        = dir('BB*');
cluster_suffix = '';
task_block     = 12;

cd(fullfile(path_manage, 'analyses_output'));

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
imgs_se_3T_0pt2          = string([]);
imgs_se_7T_0pt2          = string([]);

for folder_id = 1:length(folders)
   
   img = folders(folder_id).name;
   
   %-this analysis is only for task data
   %-for subjects 6+7 m-sequence was used
   %-one image for TR=350ms
   if ~contains(img, '_18s_off_12s_on_') || contains(img, 'TR350')
      continue
   end
   
   %-this subject was partially sleeping, at least at the end (for SE)
   if contains(img, 'BB_010_7T')
      continue
   end
   
   %-at least BB_001_3T
   if contains(img, 'MB1')
      continue
   end
   
   %-for this subject two SE BOLD tasks
   if strcmp(img, 'BB_003_3T_cmrr_mbep2d_se_task_18s_off_12s_on_20200710180104_11')
      continue
   end
   
   %-additional exclusions (artifacts)
   %{
   7:  BB 004 7T cmrr mbep2d diff b0pt2 1 task 18s off 12s on
   17: BB 004 3T cmrr mbep2d diff 0pt2 1 task 18s off 12s on
   46: BB 014 7T cmrr mbep2d diff b0pt2 1 task 18s off 12s on
   50: BB 016 7T cmrr mbep2d diff b0pt2 1 task 18s off 12s on
   %}
   if strcmp(img, 'BB_004_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200706151255_16') || ...
         strcmp(img, 'BB_004_3T_cmrr_mbep2d_diff_0pt2_1_task_18s_off_12s_on_20200804162733_15') || ...
         strcmp(img, 'BB_014_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201029092350_14_combined') || ...
         strcmp(img, 'BB_016_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20201030143609_9_combined')
      continue
   end
   
   if isempty(dir(fullfile(img, 'response_*.mat')))
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
      
      if     contains(img, '_3T_') && contains(img, 'b0_1')
         imgs_diff_0_1_3T(end+1)         = img;
      elseif contains(img, '_7T_') && contains(img, 'b0_1')
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
      
      if     contains(img, '_3T_') && contains(img, '0pt2_1') && ~contains(img, '_monopolar_')
         imgs_se_3T_0pt2(end+1)          = img;
      elseif contains(img, '_7T_') && contains(img, '0pt2_1') && ~contains(img, '_monopolar_')
         imgs_se_7T_0pt2(end+1)          = img;
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
imgs_se_3T_0pt2          = imgs_se_3T_0pt2';
imgs_se_7T_0pt2          = imgs_se_7T_0pt2';


for VOI_id = 1:4
   
   if VOI_id == 1
      VOI = 'all';
   elseif VOI_id == 2
      VOI = 'visual';
   elseif VOI_id == 3
      VOI = 'motor';
   else
      VOI = 'somatosensory';
   end
   
   figure('visible', 'off', 'PaperPositionMode', 'auto', 'PaperOrientation','landscape', 'pos', [0 0 1000 400]);
   
   for subplot_id = 1:10
      
      subplot(2, 5, subplot_id);
      
      if subplot_id == 1
         imgs    = imgs_se_3T;
         title_t = 'SE BOLD 3T';
      elseif subplot_id == 2
         imgs    = imgs_diff_0pt2_1_bi_3T;
         title_t = 'ADC 0.2/1 bipolar 3T';
      elseif subplot_id == 3
         imgs    = imgs_se_3T_0pt2;
         title_t = 'pseudo SE BOLD 3T';
      elseif subplot_id == 4
         imgs    = imgs_diff_0_1_3T;
         title_t = 'ADC 0/1 bipolar 3T';
      elseif subplot_id == 5
         imgs    = imgs_diff_0pt2_1_mono_3T;
         title_t = 'ADC 0.2/1 monopolar 3T';
      elseif subplot_id == 6
         imgs    = imgs_se_7T;
         title_t = 'SE BOLD 7T';
      elseif subplot_id == 7
         imgs    = imgs_diff_0pt2_1_bi_7T;
         title_t = 'ADC 0.2/1 bipolar 7T';
      elseif subplot_id == 8
         imgs    = imgs_se_7T_0pt2;
         title_t = 'pseudo SE BOLD 7T';
      elseif subplot_id == 9
         imgs    = imgs_diff_0_1_7T;
         title_t = 'ADC 0/1 bipolar 7T';
      elseif subplot_id == 10
         imgs    = imgs_diff_0pt2_1_mono_7T;
         title_t = 'ADC 0.2/1 monopolar 7T';
      end
      
      response_avg   = 0;
      No_of_images   = 0;
      subject_names  = cellstr({''});
      patch([0 task_block task_block 0], [-100000 -100000 100000 100000], [200 200 200]/255, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
      xlim([-8 30]);
      responses_norm = NaN(length(imgs), 38001);
      ylim_min =  Inf;
      ylim_max = -Inf;
      
      for img_id = 1:length(imgs)
         
         cd(imgs(img_id));
         run(fullfile(path_manage, 'omit_imgs_with_little_activation.m'));
         
         if (omit == false) && exist(['response_' VOI cluster_suffix '.mat'], 'file') == 2
            
            if subplot_id == 3 || subplot_id == 8
               load(['response_' VOI '_map_0pt2' cluster_suffix '.mat']);
            else
               load(['response_' VOI cluster_suffix '.mat']);
            end
            
            %-normalizing so that at the beginning of the stimulus there is value 1
            %-normalizing destroys the quantitative nature of ADC!
            response = response / response(8001);
            responses_norm(img_id,:) = response;
            
            if min(response) < ylim_min
               ylim_min = min(response);
            end
            if max(response) > ylim_max
               ylim_max = max(response);
            end
            
            time_axis = -8:0.001:30;
            plot(time_axis, response); hold on;
            response_avg = response_avg + response;
            No_of_images = No_of_images + 1;
            
            %-needed for the legend
            img = convertStringsToChars(imgs(img_id));
            img = strrep(img, '_', ' ');
            subject_names(end+1) = cellstr(img(1:9));
            
         end
         
         cd ..
         
      end
      
      %-average across subjects
      if No_of_images > 0
         response_avg = response_avg / No_of_images;
         if size(response_avg, 1) ~= 1
            response_avg = response_avg';
         end
         plot(time_axis, response_avg, 'LineWidth', 1.75, 'color', 'black');
      end
      
      if ylim_min == Inf
         ylim_min = -10;
      end
      if ylim_max == -Inf
         ylim_max = 10;
      end
      
      ylim([ylim_min ylim_max]);
      set(gca, 'FontSize', 7);
      title({title_t; [num2str(No_of_images) ' subjects']; ' '});
      xlabel('time [s]');
      
      if ~contains(title_t, 'ADC')
         ylabel('arbitrary units');
      else
         ylabel('\mum^2 / ms');
      end
      
      leg = legend(subject_names(2:end), 'FontSize', 3, 'NumColumns', 2);
      leg.ItemTokenSize = [15, 18];
      legend boxoff;
      
   end
   
   cd ../figures
   %-the option '-painters' is important as with too many lines vectorisation was vanishing
   print(['group_response_functions_' VOI], '-painters', '-dpdf');
   system(['pdfcrop group_response_functions_' VOI '.pdf group_response_functions_' VOI '.pdf']);
   cd ../analyses_output
   
end


responses_8 = NaN(8, 38001);
figure('visible', 'off', 'PaperPositionMode', 'auto', 'PaperOrientation','landscape', 'pos', [0 0 1000 800]);

for VOI_id = 1:4
   
   if VOI_id == 1
      VOI = 'all';
   elseif VOI_id == 2
      VOI = 'visual';
   elseif VOI_id == 3
      VOI = 'motor';
   else
      VOI = 'somatosensory';
   end
   
   for subplot_id = [6 1 7 2 8 3 9 4 10 5]
      
      if subplot_id < 6
         subplot(4, 5, (VOI_id-1)*5 + subplot_id);
      else
         subplot(4, 5, (VOI_id-1)*5 + subplot_id - 5);
      end
      
      if subplot_id == 1
         imgs    = imgs_se_3T;
         title_t = '(A)  SE BOLD 3T';
      elseif subplot_id == 2
         imgs    = imgs_diff_0pt2_1_bi_3T;
         title_t = '(B)  ADC 0.2/1 bipolar 3T';
      elseif subplot_id == 3
         imgs    = imgs_se_3T_0pt2;
         title_t = '(C)  pseudo SE BOLD 3T';
      elseif subplot_id == 4
         imgs    = imgs_diff_0_1_3T;
         title_t = '(D)  ADC 0/1 bipolar 3T';
      elseif subplot_id == 5
         imgs    = imgs_diff_0pt2_1_mono_3T;
         title_t = '(E)  ADC 0.2/1 monopolar 3T';
      elseif subplot_id == 6
         imgs    = imgs_se_7T;
         title_t = '(A)  SE BOLD 7T';
      elseif subplot_id == 7
         imgs    = imgs_diff_0pt2_1_bi_7T;
         title_t = '(B)  ADC 0.2/1 bipolar 7T';
      elseif subplot_id == 8
         imgs    = imgs_se_7T_0pt2;
         title_t = '(C)  pseudo SE BOLD 7T';
      elseif subplot_id == 9
         imgs    = imgs_diff_0_1_7T;
         title_t = '(D)  ADC 0/1 bipolar 7T';
      elseif subplot_id == 10
         imgs    = imgs_diff_0pt2_1_mono_7T;
         title_t = '(E)  ADC 0.2/1 monopolar 7T';
      end
      
      response_avg = 0;
      
      No_of_images = 0;
      if subplot_id > 5
         patch([0 task_block task_block 0], [-100000 -100000 100000 100000], [200 200 200]/255, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
      end
      
      xlim([-8 30]);
      responses      = NaN(length(imgs), 38001);
      responses_norm = NaN(length(imgs), 38001);
      
      for img_id = 1:length(imgs)
         
         cd(imgs(img_id));
         run(fullfile(path_manage, 'omit_imgs_with_little_activation.m'));
         
         if (omit == false) && exist(['response_' VOI cluster_suffix '.mat'], 'file') == 2
            
            if subplot_id == 3 || subplot_id == 8
               load(['response_' VOI '_map_0pt2' cluster_suffix '.mat']);
            else
               load(['response_' VOI cluster_suffix '.mat']);
            end
            if size(responses, 2) ~= length(response)
               responses = NaN(length(imgs), length(response));
            end
            responses(img_id, :) = response;
         else
            disp(['little activation or no response_' VOI '.mat for ' char(imgs(img_id))]);
         end
         
         cd ..
         
      end
      
      for img_id = 1:length(imgs)
         
         cd(imgs(img_id));
         run(fullfile(path_manage, 'omit_imgs_with_little_activation.m'));
         
         if (omit == false) && exist(['response_' VOI cluster_suffix '.mat'], 'file') == 2
            
            if subplot_id == 3 || subplot_id == 8
               load(['response_' VOI '_map_0pt2' cluster_suffix '.mat']);
            else
               load(['response_' VOI cluster_suffix '.mat']);
            end
            %-normalizing so that in the 8s window before stimulus start there is average signal of 1
            %-normalizing destroys the quantitative nature of ADC!
            response = response / mean(response(1:8000));
            responses_norm(img_id,:) = response;
            time_axis = -8:0.001:30;
            response_avg = response_avg + response;
            No_of_images = No_of_images + 1;
            
         end
         
         cd ..
         
      end
      
      %-average across subjects
      if No_of_images > 0
         response_avg = response_avg / No_of_images;
         if size(response_avg, 1) ~= 1
            response_avg = response_avg';
         end
         if contains(imgs(1), '_3T_')
            plot(time_axis, response_avg, 'LineWidth', 1.75, 'color', 'blue');
         else
            plot(time_axis, response_avg, 'LineWidth', 1.75, 'color', 'red');
         end
      end
      
      if     VOI_id == 2 && subplot_id == 1
         responses_8(1,:) = response_avg;
      elseif VOI_id == 2 && subplot_id == 6
         responses_8(2,:) = response_avg;
      elseif VOI_id == 2 && subplot_id == 2
         responses_8(3,:) = response_avg;
      elseif VOI_id == 2 && subplot_id == 7
         responses_8(4,:) = response_avg;
      elseif VOI_id == 3 && subplot_id == 1
         responses_8(5,:) = response_avg;
      elseif VOI_id == 3 && subplot_id == 6
         responses_8(6,:) = response_avg;
      elseif VOI_id == 3 && subplot_id == 2
         responses_8(7,:) = response_avg;
      elseif VOI_id == 3 && subplot_id == 7
         responses_8(8,:) = response_avg;
      end
      
      if length(imgs) > 1
         
         %-SEM = standard error of the mean; across subjects
         SEM = (std(responses_norm,0,1,'omitnan'))/sqrt(No_of_images);
         fill_x = [time_axis fliplr(time_axis)];
         if size(response_avg) ~= size(SEM)
            SEM = SEM';
         end
         
         fill_y = [(response_avg - SEM) fliplr(response_avg + SEM)];
         
         %-with too many points 'fill' can lead to problems, maybe related to graphics drivers
         %-https://github.com/altmany/export_fig/issues/161
         if contains(imgs(1), '_3T_')
            fill((fill_x(1:10:end))', (fill_y(1:10:end))', 'blue', 'FaceAlpha', 0.5, 'HandleVisibility','off');
         else
            fill((fill_x(1:10:end))', (fill_y(1:10:end))', 'red',  'FaceAlpha', 0.5, 'HandleVisibility','off');
         end
         pause(5);
         
      end
      
      % ylim([0.982 1.0238]);
      ylim([0.972 1.0242]);
      yline(1, 'HandleVisibility','off');
      set(gca, 'FontSize', 7);
      
      if strcmp(VOI, 'all')
         VOI_t = 'all brain';
      elseif strcmp(VOI, 'visual')
         VOI_t = 'visual cortex';
      elseif strcmp(VOI, 'motor')
         VOI_t = 'motor cortex';
      else
         VOI_t = 'somatosensory cortex';
      end
      
      if     subplot_id == 3 && VOI_id == 1
         title({title_t(1:end-3); ' '; VOI_t; ' '});
      elseif subplot_id ~= 3 && VOI_id == 1
         title({title_t(1:end-3); ' '; ' '; ' '});
      elseif subplot_id == 8
         title({VOI_t; ' '});
      end
      
      if subplot_id < 6
         %-be careful with what happens when 3T or 7T missing, but the other one is not missing!
         leg = legend({['7T: ' num2str(No_of_images_previous) ' subject/s'];
            ['3T: ' num2str(No_of_images)          ' subject/s']}, 'FontSize', 7, 'Location', 'Best');
         leg.ItemTokenSize = [15, 18];
         legend boxoff;
      end
      
      xlabel('time [s]');
      if ~contains(title_t, 'ADC')
         ylabel('relative BOLD signal');
      else
         ylabel('relative ADC signal');
      end
      No_of_images_previous = No_of_images;
      
   end
   
end

cd ../figures
print('group_response_functions_mean_se', '-dpdf');
print('group_response_functions_mean_se', '-dsvg');
pause(7);
system('pdfcrop group_response_functions_mean_se.pdf group_response_functions_mean_se.pdf');


figure('visible', 'off', 'PaperPositionMode', 'auto', 'PaperOrientation','landscape', 'pos', [0 0 780 250]);
for subplot_id = 1:2
   
   subplot(1, 3, subplot_id);
   patch([0 task_block task_block 0], [-100000 -100000 100000 100000], [200 200 200]/255, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
   h1 = plot(time_axis, responses_8((subplot_id-1)*4+1, :),       'color', [0,      0.5,    0     ]); hold on;
   h2 = plot(time_axis, responses_8((subplot_id-1)*4+2, :), '--', 'color', [0,      0.5,    0     ]); hold on;
   h3 = plot(time_axis, responses_8((subplot_id-1)*4+3, :),       'color', [0.4940, 0.1840, 0.5560]); hold on;
   h4 = plot(time_axis, responses_8((subplot_id-1)*4+4, :), '--', 'color', [0.4940, 0.1840, 0.5560]);
   if subplot_id == 1
      title({'visual cortex'; ''});
   else
      title({'motor cortex';  ''});
   end
   xlim([-8 30]);
   ylim([0.972 1.0242]);
   xlabel('time [s]');
   ylabel({''; 'relative BOLD/ADC signal'; ''});
   if subplot_id == 1
      leg = legend([h1 h2 h3 h4], {'SE BOLD 3T'; 'SE BOLD 7T'; 'ADC 0.2/1 bipolar 3T'; 'ADC 0.2/1 bipolar 7T'}, 'FontSize', 7, 'Location', 'Best');
      leg.ItemTokenSize = [15, 18];
      legend boxoff;
   end
   
end

print('group_response_functions_mean_8', '-dpdf');
pause(7);
system('pdfcrop group_response_functions_mean_8.pdf group_response_functions_mean_8.pdf');

cd ..
