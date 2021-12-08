function [ SampleMask ]= genSamplingMask(imgDim,pctTotal,pctCenter,pdfType, pdfParam)
% function [ SampleMask ]=
% genSamplingMask(imgDim,pctTotal,pctCenter,pdfType, pdfParam)
% 
% genSamplingMask.m Generates a variable density K-space undersampling
% mask for random Cartesian phase encoding sampling.
%
% INPUT PARAMETERS:
%
%   1) imgDim -- 2D array containing the size of the image. Example:
%   [256,256]
% 
%   2) pctTotal -- percentage of undersampleing. Example: 0.25 =>
%   samples only 25% of K-space.
% 
%   3) pctCenter -- Percentage of the sampling points that will be
%   taken from the center of K-Space.
% 
%   4) pdfType -- Type of pdf to choose k-space samples from. Choose
%   from 'Normal' (Gaussian) or Laplace.
% 
%   5) pdfParam -- distribution parameters. For Normal distribution,
%   must specify Mu (mean) and Sigma (variance = Sigma^2). For Laplace
%   distribution, must specify Mu (mean) and b (variance = 2b^2).
%   
%--------------------------------------------------------------------------
imgSizeX = imgDim(1);
kSpacePoints = 1:imgSizeX; % K-space points to sample from.
if( strcmp(pdfType,'Normal') )
    PDF = pdf('Normal', kSpacePoints, pdfParam.Mu, pdfParam.Sigma);
    PDF_norm = PDF/max(PDF(:));
elseif( strcmp(pdfType,'Laplace') )
    PDF = laplacePDF( kSpacePoints, pdfParam.Mu, pdfParam.b );
    PDF_norm = PDF/max(PDF(:));
end

numSamplePoints = round(imgSizeX*pctTotal);
numCenterPoints = floor(numSamplePoints*pctCenter);
sampleIndices = nan(numSamplePoints,1); % stores vertical k-space frequencies that will be sampled.
%% Get Center Samples
% The following code gets the center indices. It also adjusts for even/odd
% case so that the number of center points is symmeteric about the median.
if (numCenterPoints > 0)
    middleSample = median(1:imgSizeX);
    if( mod(imgSizeX,2)==0) % If dimension of K-space is even
        if(mod(numCenterPoints,2)~=0) %If number of center points is odd.
            numCenterPoints = numCenterPoints +1;
        end
        sampleIndices(1:numCenterPoints) = [ (floor(middleSample)-numCenterPoints/2 + 1):floor(middleSample); ...
                                              ceil(middleSample):(ceil(middleSample)+ numCenterPoints/2 -1)];
    else % If dimension of K-space is odd
        if(mod(numCenterPoints,2)==0) %If number of center points is even.
            numCenterPoints = numCenterPoints + 1;
        end
        sampleIndices(1:numCenterPoints) = [ (middleSample - (numCenterPoints-1)/2): (middleSample-1); middleSample; ...
                                             (middleSample+1): (middleSample + (numCenterPoints-1)/2) ];
    end 
    centerPts =  sampleIndices(1:numCenterPoints);
    PDF_norm(centerPts) = 1;
else
    centerPts = [];
end
%% Draw from PDF remaining samples.
numRemainingPoints = numSamplePoints - numCenterPoints;
if (numCenterPoints > 0)
    remainingIndices  = [1:(sampleIndices(1)-1), (sampleIndices(numCenterPoints)+1):imgSizeX]';
else
    remainingIndices = 1:imgSizeX;
end
additionalSampleIndices = datasample(remainingIndices, numRemainingPoints, 'Replace',false, 'Weights', PDF_norm(remainingIndices) );
sampleIndices(numCenterPoints+1:end) = additionalSampleIndices;
%% Make sampling mask.
% The sampling mask convention I am using has the center of K-space
% located at the corners of the mask.
mask = false(imgSizeX);
mask(sampleIndices,:) = 1;
mask = ifftshift(mask);

SampleMask.Mask = mask;
SampleMask.pctTotal = pctTotal;
SampleMask.DistributionType = pdfType;
SampleMask.PDF_norm = PDF_norm;
SampleMask.pctCenter = pctCenter;

if( strcmp(pdfType,'Normal') )
    SampleMask.DistributionParams.Mu =  pdfParam.Mu;
    SampleMask.DistributionParams.Sigma = pdfParam.Sigma;
elseif( strcmp(pdfType,'Laplace') )
    SampleMask.DistributionParams.Mu =  pdfParam.Mu;
    SampleMask.DistributionParams.b = pdfParam.b;
end

end %function

