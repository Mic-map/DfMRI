function plot_diagnostic_plots(path_to_the_map, map_name, mask_name, TR, cutoff_freq, assumed_exper_freq, true_exper_freq)

%-By Wiktor Olszowy, University of Cambridge, EPFL, OlszowyW@gmail.com
%-Written following study 'Accurate autocorrelation modeling substantially improves fMRI reliability'
%-https://www.nature.com/articles/s41467-019-09230-w.pdf
%-March 2019 - June 2020
%-This script plots diagnostic plots, e.g. of GLM residuals.
%-For GLM residuals:
%-If there is strong structure visible in the GLM residuals (the power spectra are not flat),
%-the first level results are likely confounded.

%-specify the default values for the cutoff frequency used by the high-pass filter,
%-for the assumed experimental design frequency and for the true experimental design frequency;
%-10 chosen, as it is beyond the plotted frequencies
if ~exist('cutoff_freq',        'var'), cutoff_freq        = 10; end
if ~exist('assumed_exper_freq', 'var'), assumed_exper_freq = 10; end
if ~exist('true_exper_freq',    'var'), true_exper_freq    = 10; end

%-Fast Fourier Transform (FFT) will pad the voxel-wise time series to that length with trailing zeros
%-(if no. of time points lower) or truncate to that length (if no. of time points higher)
fft_n = 512;

initial_path = pwd;
cd(path_to_the_map)

if exist([map_name '.nii.gz'], 'file') == 2
   system(['gunzip ' map_name '.nii.gz']);
end

%-read the map

if contains(map_name, 'time_course_mean')
   load(map_name);
   map = zeros([1 1 1 length(time_course_mean)]);
   map(1,1,1,:) = time_course_mean;
else
   map = niftiread([map_name '.nii']);
end

%-compressing to save space
if exist([map_name '.nii'], 'file') == 2
   system(['gzip ' map_name '.nii']);
end

%-subtract mean
map = map - mean(map, 4);

%-masking (if a mask provided)

if ~strcmp(mask_name, 'none')
   
   mask = niftiread(fullfile(initial_path, [mask_name '.nii.gz']));
   
else
   
   mask = ones(size(map));
   mask = mask(:,:,:,1);
   
end


%-calculate the power spectra

dims               = size(map);
power_spectra      = zeros(fft_n, 1);
no_of_brain_voxels = 0;

for i1       = 1:dims(1)
   
   for i2    = 1:dims(2)
      
      for i3 = 1:dims(3)
         
         ts  = squeeze(map(i1, i2, i3, :));
         
         if sum(isnan(ts)) == 0
            
            if (std(ts)    ~= 0) && (mask(i1, i2, i3) == 1)
               
               %-make signal variance equal to 1
               ts                 = ts/(std(ts) + eps);
               
               %-compute the discrete Fourier transform (DFT)
               DFT                = fft(ts, fft_n);
               
               power_spectra      = power_spectra + ((abs(DFT)).^2)/min(dims(4), fft_n);
               
               no_of_brain_voxels = no_of_brain_voxels + 1;
               
            end
            
         end
         
      end
      
   end
   
end

if no_of_brain_voxels > 0
   
   %-average power spectra over all brain voxels
   power_spectra = power_spectra / no_of_brain_voxels;
   
else
   
   disp('no signal!');
   return
   
end


figure('rend', 'painters', 'pos', [0 0 400 400], 'Visible', 'off');

%-plot the power spectra

subplot(2, 1, 1);

f      = linspace(0, 0.5/TR, 257);
max_y  = max(power_spectra);
h1     = plot(f, power_spectra(1:257), 'r'); hold on;
h4     = plot(cutoff_freq,        0,     'k*');  hold on;
h5     = plot(cutoff_freq,        max_y, 'k*');  hold on;
h6     = plot(assumed_exper_freq, 0,     'c*');  hold on;
h7     = plot(assumed_exper_freq, max_y, 'c*');  hold on;
h8     = plot(true_exper_freq,    0,     'm*');  hold on;
h9     = plot(true_exper_freq,    max_y, 'm*');  hold on;
h10    = plot([0 0.5/TR], [1 1],         'k--'); hold on;
hx     = xlabel({' ', 'Frequency [Hz]', ' '});

if contains(map_name, 'time_course_mean')
   htitle = title({map_name; ' '; 'Power spectra'}, 'interpreter', 'none');
else
   htitle = title('Power spectra', 'interpreter', 'none');
end

xlim([0 0.5/TR]);
ylim([0 max_y]);

set([h1 h10],    'LineWidth', 1.25);
set([hx htitle], 'FontSize',  7);
set(gca, 'XTick',      linspace(0, 0.5/TR, 6));
set(gca, 'XTickLabel', round(linspace(0, 0.5/TR, 6), 2));
set(gca, 'FontSize', 7);

if (power_spectra(257) < 0) || (max_y > 2)
   legend([h10 h4 h6 h8], {'Ideal power spectra', 'High pass filter frequency cutoff', 'Assumed design frequency', 'True design frequency'}, 'box', 'off', 'FontSize', 7, 'Location', 'northeast');
else
   legend([h10 h4 h6 h8], {'Ideal power spectra', 'High pass filter frequency cutoff', 'Assumed design frequency', 'True design frequency'}, 'box', 'off', 'FontSize', 7, 'Location', 'southeast');
end


%-make QQ-plot

subplot(2, 1, 2);

%-normalize the signal in each voxel separately
map = zscore(map, [], 4);
map = reshape(map, 1, []);
map = map(map~=0);
map = map(~isnan(map));
map = zscore(map);

%-control random number generation
rng('default');
rng(1);
N      = sort(normrnd(0, 1, [1 length(map)]));
h1     = plot(N, sort(map), 'color', 'r'); hold on;
h4     = plot([-1000 1000], [-1000 1000], 'k--'); hold on;
hx     = xlabel({' ', 'Normal theoretical quantiles', ' '});
hy     = ylabel('Data quantiles', 'Units', 'normalized');
htitle = title('QQ-plot', 'interpreter', 'none');
min_x  = min(N);
max_x  = max(N);
min_y  = min(map);
max_y  = max(map);
xlim([min_x max_x]);
ylim([min_y max_y]);
set([h1 h4], 'LineWidth', 1.5);
set(gca, 'FontSize', 7);
set([hx hy htitle], 'FontSize', 7);


if contains(map_name, 'time_course_mean')
   print(['diagnostic_plots_' map_name], '-dpdf');
   system(['pdfcrop diagnostic_plots_' map_name '.pdf diagnostic_plots_' map_name '.pdf']);
else
   print(['diagnostic_plots_' map_name], '-dpng');
end

cd(initial_path);


end
