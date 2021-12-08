function [ I_noisy, noise, sigma_n ] = addNoise(I, SNR)
%addNoise(I, SNR) adds Gaussian noise to an array. If the array is
%complex, complex noise is added appropriately.
%
% If the noise, n, is complex than n = x+jy where x, y are independent
% 0-mean real Gaussian each with variance sigma^2/2 so that var(n) =
% sigma^2.
%
% Input Parameters: I = original uncorrupted image.
% 
% SNR = 20*log10(sqrt( mean(mean (abs(I).^2) ) )/sigma_noise ).
%--------------------------------------------------------------------------
sigma_n = sqrt( mean( abs(I(:)).^2 )  ) / ( 10^(SNR/20) );
if ( isreal(I) )
    noise = sigma_n * randn(size(I));
else
    noise = (sigma_n/sqrt(2)) * [randn(size(I)) + j*randn(size(I)) ];
end
I_noisy = I + noise;

end %function

