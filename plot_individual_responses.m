

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Plotting mean time-courses and response functions for one exemplary subject at 3T
%%%%   and for one exemplary subject at 7T.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       April 2021 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage          = '/home/wiktor/Desktop/Dfmri/pipeline';
path_manage_software = fullfile(path_manage, 'software');
task_block           = 12;


addpath(fullfile(path_manage_software, 'extra_functions'));

cd(path_manage);

figure('visible', 'off');

for subject_id = 1:2
   
   if subject_id == 1
      cd(fullfile(path_manage, '/analyses_output/BB_013_3T_cmrr_mbep2d_diff_2pt5iso_b0pt2_1_1206_task_18s_off_12s_on_20201023162033_11_combined'));
      img_name   = 'exemplary subject at 3T: ADC 0.2/1 bipolar';
      TR         = 1;
      No_volumes = 599;
      No_act_vox = [49 28];
   elseif subject_id == 2
      cd(fullfile(path_manage, '/analyses_output/BB_003_7T_cmrr_mbep2d_diff_b0pt2_1_task_18s_off_12s_on_20200702093313_18'));
      img_name   = 'exemplary subject at 7T: ADC 0.2/1 bipolar';
      TR         = 2.07;
      No_volumes = 300;
      No_act_vox = [72 40];
   end
   
   exper_design = readtable('exper_design.txt', 'Delimiter', ' ', 'ReadVariableNames', false);
   exper_design = table2array(exper_design);
   
   subplot(2, 3, (subject_id-1)*3+1);
   load time_course_mean_visual;
   plot((1:No_volumes)*TR, time_course_mean); hold on;
   ylim_min = min(time_course_mean);
   ylim_max = max(time_course_mean);
   for i = 1:size(exper_design, 1)
      h_p = patch(exper_design(i,1)+[0 task_block task_block 0], [-1000 -1000 1000 1000], [200 200 200]/255, 'EdgeColor', 'none');
   end
   h_v = plot((1:No_volumes)*TR, time_course_mean, 'color', [0, 0.4470, 0.7410]);
   load time_course_mean_motor;
   hold on;
   h_m = plot((1:No_volumes)*TR, time_course_mean, 'color', [0.8500, 0.3250, 0.0980]);
   ylim_min = min([ylim_min min(time_course_mean)]);
   ylim_max = max([ylim_max max(time_course_mean)]);
   % legend([h_v h_m], 'visual cortex', 'motor cortex', 'AutoUpdate','off');
   % legend boxoff
   title({' '; ' '; ' '; 'time-courses'; ' '}, 'Interpreter', 'latex');
   xlim([0        No_volumes*TR]);
   ylim([ylim_min ylim_max]);
   xlabel('time [s]', 'Interpreter', 'latex');
   set(gca, 'FontSize', 5);
   
   h = subplot(2, 3, (subject_id-1)*3+2);
   load response_visual;
   plot(-8:0.001:30, response);
   patch([0 task_block task_block 0], [min(response) min(response) max(response) max(response)], [200 200 200]/255, 'EdgeColor', 'none');
   hold on;
   plot(-8:0.001:30, response, 'color', [0, 0.4470, 0.7410]);
   title({' '; img_name; ' '; 'response function: visual cortex'; ' '}, 'Interpreter', 'latex');
   xlabel('time [s]', 'Interpreter', 'latex');
   xlim([-8 30]);
   ylim([min(response) max(response)]);
   xticks([0 10 20 30]);
   xticklabels({'0','10','20','30'});
   annotation('textbox', 'String', [num2str(No_act_vox(1)) ' active voxels'], 'Position', h.Position, 'Vert', 'bottom', 'Horizontal', 'left', 'EdgeColor', 'none', 'FontSize', 5, 'Interpreter', 'latex');
   set(gca, 'FontSize', 5);
   
   h = subplot(2, 3, (subject_id-1)*3+3);
   load response_motor;
   plot(-8:0.001:30, response, 'color', [0.8500, 0.3250, 0.0980]);
   patch([0 task_block task_block 0], [min(response) min(response) max(response) max(response)], [200 200 200]/255, 'EdgeColor', 'none');
   hold on;
   plot(-8:0.001:30, response, 'color', [0.8500, 0.3250, 0.0980]);
   title({' '; ' '; ' '; 'response function: motor cortex'; ' '}, 'Interpreter', 'latex');
   xlabel('time [s]', 'Interpreter', 'latex', 'Interpreter', 'latex');
   xlim([-8 30]);
   ylim([min(response) max(response)]);
   xticks([0 10 20 30]);
   xticklabels({'0','10','20','30'});
   annotation('textbox', 'String', [num2str(No_act_vox(2)) ' active voxels'], 'Position', h.Position, 'Vert', 'bottom', 'Horizontal', 'left', 'EdgeColor', 'none', 'FontSize', 5, 'Interpreter', 'latex');
   set(gca, 'FontSize', 5);
   
end

cd(fullfile(path_manage, 'figures'));
print('individual_responses', '-dpdf');
pause(7);
system('pdfcrop individual_responses.pdf individual_responses.pdf');
