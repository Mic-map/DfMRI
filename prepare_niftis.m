

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   This code needs to be run with caution. After data acquisition, convert DICOMs to NIfTIs,
%%%%   ideally with MRIcron, and look at the NIfTIs, e.g. with fsleyes. Be sure you understand
%%%%   how many volumes were acquired for which b-value, and in which order. This code is used
%%%%   to do initial preprocessing of some of the datasets.
%%%%
%%%%   There is quite a lot of variability across the datasets in terms of the exact order of
%%%%   b-values. Some acqusitions are interleaved (0.2/1/0.2/1/0.2/1/...), while some other
%%%%   acquisitions are primarily for one b-value. (In the latter case, the ADC maps are later
%%%%   calculated based on two acquisition runs.) However, if the acquisition is primarily for
%%%%   one b-value, there are still some volumes for b-value 0. Such volumes are acquired,
%%%%   e.g., to be able to assess drifts in the data or to enable better calibration of the
%%%%   system.
%%%%
%%%%   This code:
%%%%     (1) Selects revPE volumes with the lower b-value (in almost all cases = 0.2), as
%%%%         these images have higher SNR, so better topup later. The output has suffix '_short'.
%%%%     (2) Interpolates signal in the middle volume, where b-value = 0. This middle b-value 0
%%%%         volume was only acquired in some sessions (always first check the NIfTI with fsleyes!).
%%%%     (3) Removes the last volume, where b-value = 0. Again, this is needed only in some
%%%%         situations (always first check the NIfTI with fsleyes!).
%%%%     (4) Merges the volumes from the two separate acquisition runs into one NIfTI: suffix
%%%%         '_combined'.
%%%%
%%%%   Please mind there are differences between 3T and 7T (comment/uncomment).
%%%%
%%%%   'img*_name' variables have to be initiated manually.
%%%%
%%%%   After you run this code on the images you wanted to process, visually check all the new
%%%%   NIfTIs, e.g. in the movie-mode in fsleyes. If some mistake is made here, the later
%%%%   analysis of the data is very problematic.
%%%%
%%%%   Remember that in fsleyes numbering starts with 0! For movie-mode, turn off 'Synchronise
%%%%   movie updates' (settings behind the "monkey wrench" button).
%%%%
%%%%   For each subject, copy all the relevant NIfTIs, including MP(2)RAGE and other datasets
%%%%   which did not need to be preprocessed with this code, to '/pipeline/scans/'.
%%%%
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       October 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd /home/wiktor/Desktop/Dfmri/scans_all/scans_all_niftis/2021_02_26_BB_024_7T/


for img_id = 1:3
   
   if img_id == 1
      %-needs to be changed!!!!!!!
      img_name = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_rs_revPE_monopolar_20210226150433_9';
   elseif img_id == 2
      %-needs to be changed!!!!!!!
      img_name = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_monopolar_20210226150433_7';
   elseif img_id == 3
      %-needs to be changed!!!!!!!
      img_name = 'BB_024_7T_cmrr_mbep2d_diff_b0pt2_task_revPE_bipolar_20210226150433_8';
   end
   
   img_hdr = niftiinfo(img_name);
   img = niftiread(img_hdr);
   %-difference between 7T and 3T results from different vectors of diffusion encoding used at the two systems (for revPE) - different scaling -> different ordering of b-values
   %-for 7T
   img = img(:,:,:,4:end);
   %-for 3T
   %	img = img(:,:,:,3:2:end);
   img_hdr.ImageSize(4) = size(img, 4);
   niftiwrite(img, [img_name '_short'], img_hdr, 'compressed', 1);
   
end


for img_id = 1
   
   if img_id == 1
      %-needs to be changed!!!!!!!
      img1_name = 'BB_021_7T_cmrr_mbep2d_diff_b0pt2_task_18s_off_12s_on_monopolar_20210219110948_10';
      img2_name = 'BB_021_7T_cmrr_mbep2d_diff_1_task_18s_off_12s_on_monopolar_20210219110948_11';
   end
   
   %-I am changing the format, as the int format led to some scaling problems
   system(['fslmaths ' img1_name ' ' img1_name ' -odt float']);
   system(['fslmaths ' img2_name ' ' img2_name ' -odt float']);
   
   img1_hdr  = niftiinfo(img1_name);
   img2_hdr  = niftiinfo(img2_name);
   img1      = niftiread(img1_name);
   img2      = niftiread(img2_name);
   
   %-difference between 7T and 3T results from Prisma enabling mosaic only if additionally one b-value acquired first (before the vector of diffusion encoding is read...)
   %-for 7T
   img1(:,:,:,301) = (img1(:,:,:,300) + img1(:,:,:,302))/2;
   img2(:,:,:,301) = (img2(:,:,:,300) + img2(:,:,:,302))/2;
   
   %-for 3T
   %	img1(:,:,:,302) = (img1(:,:,:,301) + img1(:,:,:,303))/2;
   %	img2(:,:,:,302) = (img2(:,:,:,301) + img2(:,:,:,303))/2;
   hdr_combined = img1_hdr;
   hdr_combined.ImageSize(4)       = 2*hdr_combined.ImageSize(4)-2;
   
   %-for such separate acquisitions (lower b-value in one run, higher b-value in another run), the TR of the new, combined NIfTI is half the initial TR
   hdr_combined.PixelDimensions(4) = hdr_combined.PixelDimensions(4)/2;
   hdr_combined.Datatype           = 'double';
   img_combined = zeros(hdr_combined.ImageSize);
   img_combined(:,:,:,1:2:end) = img1(:,:,:,1:(end-1));
   img_combined(:,:,:,2:2:end) = img2(:,:,:,1:(end-1));
   niftiwrite(img_combined, [strrep(img1_name, '0pt2', '0pt2_1') '_combined'], hdr_combined, 'compressed', 1);
   
end
