%%
%% comb(x, k, weights)
%%
%% Apply a comb filter to some data.
%%
%% x: The data to be filtered
%% k: lag
%% weights: coeffients of x, x-k, x-2k etc...
function c = comb(x, k, weights)
	N = length(weights);
	c = zeros(1, length(x));

	for i=1:N
		c((N-1)*k+1:end)= c((N-1)*k+1:end) + weights(i) * x((N-i)*k+1 : end -(i-1)*k);
	end

