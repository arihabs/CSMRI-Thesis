function [metrics] = errorMetrics(A, B, varargin)
% function [metrics] = errorMetrics(A, B, varargin)
% 
% errorMetrics.m Computes error metrics for image comparison. 

% Compute the mean squared error (MSE), RMSE, SNR, and PSNR [1] as
% well as the MSSIM [2] of an image approximation. 

% Input Parameters:
% Image A is the reference image and image B is the approximation.
%
% Output Parameters:
% metrics -- struct containing MSE, RMSE, SNR, PSNR and MSSIM.
% 
% References:
% 
% [1] "Data Compression: The complete reference", 4th Ed.
% 
% [2] Z. Wang, A. C. Bovik, H. R. Sheikh and E. P. Simoncelli, "Image
% quality assessment: From error visibility to structural similarity,"
% IEEE Transactions on Image Processing, vol. 13, no. 4, pp. 600-612,
% Apr. 2004.
% 
% 
% Author: Ariel Habshush
% The Cooper Union for the Advancement of Science and Art,
% Department of Electrical Engineering
%
% Email: habshu@cooper.edu
% August 2013; Last revision: 18-July-2014
%--------------------------------------------------------------------------

% Check that both images are real valued.
if( any(imag(A(:))~=0) || any(imag(B(:))~=0) )
    error('Error metrics function not defined for complex valued images.')
end

if(nargin < 4)
    ParamsSSIM = [];
else
    ParamsSSIM = varargin{2};
end

if(nargin < 3)
    PixelRange = [];
else
    PixelRange = varargin{1};
end
    

% If specified, check that Pixel Range is valid (distinct numbers >= 0).
if( ~isempty(PixelRange) )
    if( any(PixelRange < 0) )
        error('Not a valid pixel range.')
    elseif( PixelRange(1) == PixelRange(2) )
        error('Not a valid pixel range.')
    elseif( PixelRange(1) > PixelRange(2) )
        PixelRange = [PixelRange(2), PixelRange(1)];
    end
    
    % Convert images back to there original dynamic range.
    A = double2int(A, PixelRange);
    B = double2int(B, PixelRange);
else
    PixelRange = [0,1];
end

maxVal = PixelRange(2);
MSE = mean( (A(:)-B(:) ).^2);
RMSE = sqrt(MSE);
PSNR = 20*log10(maxVal/RMSE );
SNR = 20*log10( sqrt( mean(A(:).^2))/RMSE );

if( ~isempty(ParamsSSIM) )
    [MSSIM, ~] = ssim_wang(A, B, ParamsSSIM.K, ParamsSSIM.window, maxVal);
else
    [MSSIM, ~] = ssim_wang(A, B, 5, 5, maxVal);
end

metrics.MSE = MSE;
metrics.RMSE = RMSE;
metrics.PSNR = PSNR;
metrics.SNR = SNR;
metrics.MSSIM = MSSIM;

end %function

