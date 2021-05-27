

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   The other codes produce figures in pdf. They are saved in the 'figures' folder. This
%%%%   code produces pngs from these pdfs.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       April 2021 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd /home/wiktor/Desktop/Dfmri/pipeline/figures/

pdfs = dir('*.pdf');

for pdf_id = 1:length(pdfs)
   pdf = pdfs(pdf_id).name;
   system(['pdftoppm ' pdf ' ' pdf(1:end-4) ' -png -r 600']);
end
