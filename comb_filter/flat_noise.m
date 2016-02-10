%%
%% flat_noise(half_length, nonflatness, F)
%%
%% Make spectrally flat noise, with a random non-flatness.
%% The magnitude of the noise in frequency space is 1. With 
%% nonflatness set to 0, the abs(fft(noise)) == [1 1 ...]. While
%% random noise is on average flat, it is of course noisy, so it can
%% be hard to pick out the shape when it's filtered.
%%
%% Parameters
%%
%% half_length: total length is 1+ 2 * half_length
%% nonflatness: noise to add to the frequency spectrum
%% F: frequency of the data
%%
%% Returns:
%% noise: spectrally flat noise
%% noise_fft_ind: the useful part of the fft of the noise
%% noise_fft_freq: frequency components of the fft of the noise
function [noise, noise_fft_ind, noise_fft_freq] = flat_noise(noise_pts, nf, F)

	fdata = exp(randn(1, noise_pts) * 2 * pi * j) .* (1 + randn(1, noise_pts)*nf);
	fdata = [1 fdata conj(fdata(end:-1:1))];
	noise = real(ifft(fdata));
	noise_fft_ind = [1:noise_pts];

	noise_fft_freq = ([0:length(noise)-1]/length(noise)) * F;
	noise_fft_freq = noise_fft_freq(noise_fft_ind);

