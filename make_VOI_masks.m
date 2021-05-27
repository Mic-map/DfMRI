

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   Making VOI (= ROI) masks. Using NMM and JHU WM atlases. Needs to be run only once -
%%%%   before the single-subject analyses.
%%%%
%%%%   Written by:    Wiktor Olszowy, CIBM Center for Biomedical Imaging, EPFL
%%%%   Contact:       olszowyw@gmail.com
%%%%   Created:       April 2020 - May 2021
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd /home/wiktor/Desktop/Dfmri/pipeline/software/brain_parcellation/


%-visual cortex
%{
<label><index>108</index><name>Right Calc calcarine cortex</name></label>
<label><index>109</index><name>Left Calc calcarine cortex</name></label>
<label><index>156</index><name>Right OCP occipital pole</name></label>
<label><index>157</index><name>Left OCP occipital pole</name></label>
<label><index>196</index><name>Right SOG superior occipital gyrus</name></label>
<label><index>197</index><name>Left SOG superior occipital gyrus</name></label>
<label><index>128</index><name>Right IOG inferior occipital gyrus</name></label>
<label><index>129</index><name>Left IOG inferior occipital gyrus</name></label>
<label><index>114</index><name>Right Cun cuneus</name></label>
<label><index>115</index><name>Left Cun cuneus</name></label>
<label><index>134</index><name>Right LiG lingual gyrus</name></label>
<label><index>135</index><name>Left LiG lingual gyrus</name></label>
<label><index>144</index><name>Right MOG middle occipital gyrus</name></label>
<label><index>145</index><name>Left MOG middle occipital gyrus</name></label>
<label><index>160</index><name>Right OFuG occipital fusiform gyrus</name></label>
<label><index>161</index><name>Left OFuG occipital fusiform gyrus</name></label>
%}

for ROI = [108 109 156 157 196 197 128 129 114 115 134 135 144 145 160 161]
   system(['fslmaths rlabels_Neuromorphometrics -thr  ' num2str(ROI-0.5) ' tmp']);
   system(['fslmaths tmp                        -uthr ' num2str(ROI+0.5) ' tmp']);
   if exist('mask_VOI_visual.nii.gz', 'file') ~= 2
      copyfile('tmp.nii.gz', 'mask_VOI_visual.nii.gz');
   else
      system('fslmaths mask_VOI_visual -add tmp mask_VOI_visual');
   end
end


%-motor cortex
%{
<label><index>150</index><name>Right MPrG precentral gyrus medial segment</name></label>
<label><index>151</index><name>Left MPrG precentral gyrus medial segment</name></label>
<label><index>182</index><name>Right PrG precentral gyrus</name></label>
<label><index>183</index><name>Left PrG precentral gyrus</name></label>
<label><index>192</index><name>Right SMC supplementary motor cortex</name></label>
<label><index>193</index><name>Left SMC supplementary motor cortex</name></label>
%}

for ROI = [150 151 182 183 192 193]
   system(['fslmaths rlabels_Neuromorphometrics -thr  ' num2str(ROI-0.5) ' tmp']);
   system(['fslmaths tmp                        -uthr ' num2str(ROI+0.5) ' tmp']);
   if exist('mask_VOI_motor.nii.gz', 'file') ~= 2
      copyfile('tmp.nii.gz', 'mask_VOI_motor.nii.gz');
   else
      system('fslmaths mask_VOI_motor -add tmp mask_VOI_motor');
   end
end


%-somatosensory cortex
%{
<label><index>148</index><name>Right MPoG postcentral gyrus medial segment</name></label>
<label><index>149</index><name>Left MPoG postcentral gyrus medial segment</name></label>
<label><index>176</index><name>Right PoG postcentral gyrus</name></label>
<label><index>177</index><name>Left PoG postcentral gyrus</name></label>
%}

for ROI = [148 149 176 177]
   system(['fslmaths rlabels_Neuromorphometrics -thr  ' num2str(ROI-0.5) ' tmp']);
   system(['fslmaths tmp                        -uthr ' num2str(ROI+0.5) ' tmp']);
   if exist('mask_VOI_somatosensory.nii.gz', 'file') ~= 2
      copyfile('tmp.nii.gz', 'mask_VOI_somatosensory.nii.gz');
   else
      system('fslmaths mask_VOI_somatosensory -add tmp mask_VOI_somatosensory');
   end
end

system('gunzip mask_VOI_visual.nii.gz');
system('gunzip mask_VOI_motor.nii.gz');
system('gunzip mask_VOI_somatosensory.nii.gz');
copyfile('rlabels_Neuromorphometrics.nii', 'mask_VOI_all.nii');
delete('tmp.nii.gz');


%-extending the NMM atlas by labels from the JHU WM atlas

system('fslmaths JHU-ICBM-labels-2mm            -add 1000 JHU-ICBM-labels-2mm_added_1000');
system('fslmaths JHU-ICBM-labels-2mm_added_1000 -thr 1001 JHU-ICBM-labels-2mm_added_1000');
system('fslmaths JHU-ICBM-labels-2mm_added_1000 -max rlabels_Neuromorphometrics NMM_JHU_atlases');
