clear

%Set up graph formatting
LW='LineWidth';
lw=2;
FS='FontSize';
Fs=12;
FN='FontName';
fn = 'Helvetica';
XLabel = @(x) xlabel(x, FS, Fs, FN, fn);
YLabel = @(x) ylabel(x, FS, Fs, FN, fn);
SetFont = @(x) set(gca, FN, fn, FS, Fs);
SetLeg =@(x) set(legend, FS, Fs, FN, fn);


%Sampling frequency of our system.
F = 300; % Hz

fs = [0:0.1:F/2];
T = 1/F;
k = 18; %This gives a 50Hz comb. 



noise_pts = 100000;
[noise, noise_fft_ind, noise_fft_freq ] = flat_noise(noise_pts, 1/50, F);

PSD = @(x) abs(fft(x));

%Theorerical comb filter response fiven some weights
comb_H = @(f, k, weights)  (exp(2*pi*T*k*j*f'*[0:length(weights)-1]) * weights')';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot of basic comb filter
%

%Weights for the basic comb filter.
basic_comb=[1 -1];


%Apply the filter to some actual data, and pad with zeros
psd = PSD(comb(noise, k, basic_comb));

clf
plot(noise_fft_freq, psd(noise_fft_ind), 'r')
hold on
plot(fs, abs(comb_H(fs, k, basic_comb)), LW, lw);



XLabel('Frequency');
YLabel('Magnitude response');
SetFont();
p = get(gcf, 'Position');
set(gcf, 'Position', [ p(1), p(2), 640, 200])
axis([0 F/2, 0, 4])
legend('Filtered signal', 'Predicted H(z)')
SetLeg();

drawnow('postscript', 'comb-1.ps');
drawnow()
p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot of next comb filter
%

%Weights for the basic comb filter.
next_comb=[1 -.5 -.5];


%Apply the filter to some actual data, and pad with zeros
psd = PSD(comb(noise, k, next_comb));

clf
plot(fs, abs(comb_H(fs, k, basic_comb)), LW, lw, 'g');
hold on
plot(noise_fft_freq, psd(noise_fft_ind), 'r')
plot(fs, abs(comb_H(fs, k, next_comb)), LW, lw);



XLabel('Frequency');
YLabel('Magnitude response');
SetFont();
p = get(gcf, 'Position');
set(gcf, 'Position', [ p(1), p(2), 640, 200])
axis([0 F/2, 0, 4])
legend('Original filter', 'Filtered signal', 'Predicted H(z)')
SetLeg();

drawnow('postscript', 'comb-2.ps');
drawnow()


optimize_weights;
func = @(x, weights)  (abs(exp(j*x'*[0:length(weights)-1]) * weights'))';

% Now draw some graphs of the filter and its phase.
clf

ol=[0:0.001:6*pi];
plot(ol, func(ol, wnew), 'g', LW, lw)
hold on
plot(omega, target, 'r', LW, lw);

legend('Filter', 'Design')
SetLeg();

axis([0, 2*pi, 0, 2.1]);
set(gca, 'XTick', [0:4]/2 * pi)
set(gca, 'XTickLabel', {'0',  '\pi/2', '\pi', '3\pi/2', '2\pi'});

%Worst case in the notch region
r = max(func(omega(find(omega<=start)), wnew));
%Centre of the notch
t = func(0, wnew);

db = @(x) 20 * log(x) / log(10)


title(sprintf('Order %i filter response. Notch depth @ Centre:  %.1f dB  (%.3e), Region (worst): %.1f dB (%.3e)', n-1, db(t), t, db(r), r), FN, fn)
XLabel('Normalised frequency (radians)')
YLabel('Magnitude')
SetFont();

drawnow('postscript', 'comb-3.ps');
drawnow()



%Apply the filter to some actual data, and pad with zeros

[noise, noise_fft_ind, noise_fft_freq ] = flat_noise(noise_pts, 0, F);

psd = PSD(comb(noise, k, wnew));

clf
plot(fs, abs(comb_H(fs, k, basic_comb)), LW, lw, 'g');
hold on
plot(noise_fft_freq, psd(noise_fft_ind), 'r')
plot(fs, abs(comb_H(fs, k, wnew)), LW, lw);



XLabel('Frequency');
YLabel('Magnitude response');
SetFont();
p = get(gcf, 'Position');
set(gcf, 'Position', [ p(1), p(2), 640, 200])
axis([0 F/2, 0, 4])
legend('1st order filter', 'Filtered signal', 'Predicted H(z)')
SetLeg();

drawnow('postscript', 'comb-4.ps');
drawnow()


for i = { '1', '2', '3', '4' }
	system([ 'pstopnm -landscape -stdout -xborder=0 -yborder=0 -xsize 4096 comb-' i{1} '.ps | pnmscale -r 4 | pnmflip -r180 | pnmtopng > comb-' i{1}' .png' ] )
end

