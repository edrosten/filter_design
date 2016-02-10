%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comb filter designer.
% 

%The filter is designed by providing a target function on [0, pi). 
%Since the filter has only cosines, this is reflected around 0, then
%repeated ad-nauseum
omega=[0:0.01:pi];

%Let the target function start at zero, then step up to 1, but with the
%step being a fragment of some function to make things smoother. Sharper transistions
%lead to larger Gibbs overshoots
%
%Note notches with a slow start (e.g. raised cosine, start > 0) require high order terms
%as features the width of the start region are required, and that's narrow. Ones that ramp
%linearly are very easy as it's very easy to force the functio to be 0 at 0.
width = 0.3;
start = 0.1;
target = ones(size(omega));
ind = find( (omega < width + start) & (omega > start));
target(ind) = .5 - .5 * cos((omega(ind)-start) / width * pi);    %Raised cosine section
%target(ind) = (omega(ind)-start) / width;                         %Straight line
target(find(omega<=start))=0;

% Set the relative importance for different regions. This allows you to trade
% accurace in one region off against another. In this case heavily weight the zero
% part of the the function (i.e. the notch)
importance = ones(size(omega));
importance(find(omega<=start))=10000000;


%Set the number of terms (this is 1 + the order)
n = 40;

% These are the absolute value of the transfer function for
% a filter with the Z transform of:
%       ___
%       \        -i k
% H(z) = >   w  z
%       /__   i
%        i
%
% with a lag of k.

func = @(x, weights)  (abs(exp(j*x'*[0:length(weights)-1]) * weights'))';
phase = @(x, weights)  (arg(exp(j*x'*[0:length(weights)-1]) * weights'))';

%Starting parameters. 
%
%Note, I'm interested in designing filters which have the coefficient of
%z positive and all powers have negative coefficients. This is required for
%filters which have f(0) = 0, and importantly, doubles the frequency of a given
%lag, meaning you double the lag size for the same frequency. This gives more
%flexibility with a low sample rate at the penalty of filtering out DC.
%
%So, choose some starting weights that (emperically) do very roughkly the 
%right thing, and scale them so the initial function is about the same scale
%as the target function.
weights = [1 -1./ 2.^[1:n-1] ];   %Decay of weights
weights(2:end) = -weights(2:end) / sum(weights(2:end)); %Notch at zero
weights = weights / max(func(omega, weights)') * max(target); % Scale


%Choose an error function which minimizes the sum of some function of the 
%absolute deviation, but weighted by the relative importance.
err = @(w) sum((abs((func(omega, w) - target)) .^6).*importance );

%Do the optimization.
wnew = fminunc(err, weights);




% Now draw some graphs of the filter and its phase.
%clf
%subplot(1, 1, 1)
%ol=[0:0.001:6*pi];
%hold on
%plot(ol, func(ol, [.5, -.5]), 'k')
%plot(ol, func(ol, weights))
%plot(ol, func(ol, wnew), 'g')
%plot(omega, target, 'r');
%
%legend('1st order', 'Initial weights', 'Optimized', 'Target');
%
%%Worst case in the notch region
%r = max(func(omega(find(omega<=start)), wnew));
%%Centre of the notch
%t = func(0, wnew);
%
%db = @(x) 20 * log(x) / log(10)
%
%
%title(sprintf('Order %i filter response. Notch depth @ Centre:  %.1f dB  (%.3e), Region (worst): %.1f (%.3e)', n-1, db(t), t, db(r), r))
%xlabel('Normalised frequency (radians)')
%ylabel('Magnitude')
%
%
%subplot(3, 1, 2)
%hold on
%
%plot(ol, phase(ol, [.5, -.5])*180/pi, 'k')
%plot(ol, phase(ol, weights) * 180 / pi)
%plot(ol, phase(ol, wnew) * 180 / pi, 'g')
%
%legend('1st order', 'Initial weights', 'Optimized');
%xlabel('Normalised frequency (radians)')
%ylabel('Phase (degrees)')
%
%
%
%subplot(3, 1, 3)
%o = [-5:0.01:1];
%semilogx(10.^o, 20 *log(func(10.^o, wnew)) / log(10))
%ylabel('Magnitude (dB)')
%xlabel('Normalised frequency')
%title('Response as a high pass filter near DC')
