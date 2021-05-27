

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Merging smaller figures into bigger ones. For each task dataset, making one pdf from
%%%%   the carpet plot, spatial zfstat1 and time courses/response functions.
%%%%   Also, merging figures across subjects, producing long pdfs - with many pages, where
%%%%   each page refers to one subject/dataset.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       October 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


path_manage = '/home/wiktor/Desktop/Dfmri/pipeline';

cd(fullfile(path_manage, 'analyses_output'));


files = dir('*/time_courses_and_response_functions_*.pdf');

for file_id = 1:length(files)
   
   disp(file_id);
   cd(files(file_id).folder);
   file       = fullfile(files(file_id).folder, files(file_id).name);
   sp_dist    = fullfile(files(file_id).folder, 'adc.feat',                           'rendered_thresh_zfstat1.png');
   
   if exist(sp_dist, 'file') ~= 2
      sp_dist = fullfile(files(file_id).folder, 'ardfMRI_denoised_unringed_sdc.feat', 'rendered_thresh_zfstat1.png');
   end
   
   if ~isempty(dir(fullfile(pwd, '*_with_sp_dist_*.pdf')))
      continue
   end
   
   carpet_plot = dir('dfMRI_raw_carpet_plot_*.pdf');
   
   system(['pdfjam ' carpet_plot.name ' ' sp_dist ' --nup 1x2 --outfile tmp.pdf']);
   system( 'pdfcrop tmp.pdf tmp.pdf');
   system(['pdfjam tmp.pdf ' file ' --nup 2x1 --landscape --outfile ' file(1:end-4) '_with_sp_dist.pdf']);
   system(['pdfcrop ' file(1:end-4) '_with_sp_dist.pdf ' file(1:end-4) '_with_sp_dist.pdf']);
   
end


cd(fullfile(path_manage, 'analyses_output'));

fig_types = cellstr({'ROI_correlations'; 'dfMRI_raw_carpet_plot'; 'time_courses_and_response_functions_*with_sp_dist'});

for fig_type_id = 1:length(fig_types)
   
   mkdir(fullfile(path_manage, '/figures/', char(fig_types(fig_type_id))));
   
   figs = dir(['*/' char(fig_types(fig_type_id)) '*']);
   
   for fig_id = 1:length(figs)
      copyfile(fullfile(figs(fig_id).folder, figs(fig_id).name), fullfile(path_manage, '/figures/', char(fig_types(fig_type_id)), '/', figs(fig_id).name));
   end
   
end
