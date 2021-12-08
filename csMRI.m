function [Images, metrics, csMRIout] = csMRI(img,kSpaceMask,Params)
% [Images, metrics, csMRIout] = csMRI(img,kSpaceMask,Params)

% Performs compressed sensing (CS) recovery techniques in order to
% reconstruct an MRI image from its undersampled k-space samples.
%
% INPUT PARAMETERS: Unless specified otherwise, all parameters listed
% are required.
% 
%       1) img -- real valued gray-scaled image to perform CS on.
% 
%       2) kSpaceMask -- binary mask indicating which points in
%       k-space will be sampled. Locations with values of 0 will not
%       be sampled. Locations with values of 1 will be sampled. The
%       mask must be the same size as img. NOTE: THE MASK MUST BE
%       ARRANGED SUCH THAT THE CENTER OF K-SPACE (CONTAINING THE MOST
%       ENERGY) IS LOCATED AT THE MASK CORNERS. IN OTHER WORDS, THE
%       MASK IS ASSUMING THAT THE K-SPACE IS GENERATED VIA calling
%       kspace = fft2(img) WITHOUT CALLING fftshift(kspace)
%       afterwards.
% 
%       3) Params -- struct whose fields contain required and optional
%       parameters for the CS-MRI algorithms.These parameters are as
%       follows:
% 
%              A) PixelRange -- a 2-element vector containing minimum
%              and maximum value a pixel in img can take. This value
%              is used to convert the pixels to double precision in
%              the range [0.0 - 1.0] for processing. The image is then
%              converted back to original dynamic range, PixelRange.
% 
%              B) numIterCsMRI -- Number of iterations to run the
%              CS-MRI algorithm. Must be an integer number >= 1.
% 
%              C) lambda -- positive constant such that nu =
%              lambda/sigma_n, where nu is the data consistency term
%              in the optimization problem and sigma_n is the standard
%              deviation of additive noise to the measurement process.
% 
%              D) PatchDim -- 2-element vector containing the
%              dimension of the patches the image is segmented into.
% 
%              E) slidingFactor -- distance between adjacent patches.
%              Maximal patch overlap occurs when slidingFactor = 1.
% 
%              F) noiseSNR -- SNR of image after white Gaussian noise
%              is added to the process. If no noise is to be added set
%              noiseSNR = Inf.
% 
%              G) ompErrorThr -- Threshold value as a stopping
%              criterion when performing orthogonal pursuit algorithm
%              (OMP) to find a sparse representation of a signal. The
%              threhold value is the Mean Squared Error of the actual
%              signal and its sparse representation.
% 
%              H) maxNZ -- Maximum number of nonzero values allowed in
%              the sparse representation of a patch.
% 
%              I) InitializationMethod -- How to initialize the
%              "sparsifying" dictionary used to represent the patch
%              signals sparsely. Choose from the following methods:
%                        a) DCT -- Use an overcomplete discrete cosine
%                        dictionary (DCT). (This function will
%                        generate the dictionary).
% 
%                        b) PCA -- Perform principle component
%                        analysis (PCA) on a subset of the image
%                        patches to find a sparse representation of
%                        the images. (This function will compute the
%                        dictionary).
% 
%                        c) External -- Provide an external
%                        dictionary. Requires the user to specify the
%                        field "initialDict" in the Params struct.
% 
%              J) InitialDict (required if InitializationMethod =
%              'External') -- contains an initialized sparsifying
%              dictionary.
% 
%              K) numAtoms (required if InitializationMethod = 'DCT'
%              or 'PCA')-- number of atoms (i.e. columns) to be used
%              in the sparsifying dictionary.
% 
%=============== K-SVD Parameters ==============
%              L) adaptDict -- When set to 1, the K-SVD algorithm will
%              be used to update the sparsifying dictionary. When set
%              to 0, sparse coding will only be performed with the
%              initialized dictionary and the dictionary will remain
%              constant througput the algorithm.
% 
%              M) numiterateKSVD (required if adaptDict = 1) -- Number
%              of iterations to run the K-SVD algorithm each time the
%              algorithm is called.
%
%==== SSIM Parameters (refer to [4] for expanded parameter
%descriptions.)
%              N) SSIM_K -- stability parameter. O) SSIM_window --
%              smoothing/averaging parameter.
%--------------------------------------------------------------------------
% 
% OUTPUT PARAMETERS:
% 
%       1) Images -- struct containing the following images:
%       
%              A) Recon -- The reconstructed image using CS-MRI
%              technique.
% 
%              B) ZF -- Zero-filled Fourier reconstructed (magnitude)
%              image. This image results when the missing K-space
%              points are filled with zeros and the signal is
%              transformed back to the image domain via an inverse FFT
%              on the zero-filled k-space data. The result is a
%              complex valued image. ZF is the magnitude of this
%              image. ZF is used as the initializtion image for the
%              CS_MRI algorithm.
% 
%              C)Noisy -- The original image + additive noise (when
%              noiseSNR == Inf, this image is not returned).
% 
%       2) Metrics -- Struct containing error metrics for the
%       following images: Recon, ZF, Noisy. For the the reconstructed
%       image, the metrics are saved after each iteration of csMRI.
% 
%       3) csMRIout -- Struct containing remaining miscellaneous
%       variables:
% 
%             A) DictInitial -- Initial sparsifying dictionary.
% 
%             B) DictFinal -- Fianl sparsifying dictionary if
%             dictionary adaptation is applied.
% 
%             C) NumNzInPatch -- Matrix contianing the number of
%             nonzero elements in each patch's sparse representation
%             for each iteration of CS-MRI. Rows correspond to
%             iteration # and columns correspond to patch #.
% 
%             D) AvgNzPerIter -- Average nonzeros elements per patch
%             per iteration. This is the mean of NumNzInPatch across
%             its columns.
% 
%             E) sigma_n -- equaivalent variance of the added Guassian
%             noise of the specified SNR.
% 
%             F) TimePerIteration -- time (in seconds) of each
%             iteration of csMRI.
% 
%             G) AvgTimePerIter -- aveage time per iteration (mean
%             value of TimePerIteration).
%==========================================================================
% ALGORITHM DESCRIPTION: Compressed Sensing MRI Via Sparse and
% Redundant Representations algorithm solves an optimization problem
% by breaking it up into stages:
%     * Break up the estimated image into patches.
% 
%     * (Optional stage) Perform K-SVD to learn and adapt the
%     sparsifying dictionary to the estimated image.
% 
%     * Perform sparse coding to determine the sparse representation
%     of each patch in the sparsifying dictionary.
% 
%     * Update the estimated image by using the K-Space update formula
%     found in DLMRI paper.
%==========================================================================
%
% Other m-files required:
%   addNoise.m, my_im2col.m (taken from [1]), overcompleteDCTdict.m,
%   pcaDictionary.m, myKSVD.m, OMPerrn.m
%==========================================================================
% REFERENCES:
% 
%   [1] S.Ravishankarand Y. Bresler, "MR image reconstruction from
%   highly undersampled k-space data by dictionary learning," IEEE
%   Trans. Med. Imag., vol. 30, no. 5, pp. 1028 - 1041, 2011.
% 
%   [2]	M. Aharon, M. Elad, and A. Bruckstein, "K-SVD: An algorithm
%   for designing overcomplete dictionaries for sparse
%   representation," IEEE Transactions on signal processing, vol. 54,
%   no. 11, pp. 4311 - 4322, 2006.
% 
%   [3] M. Aharon and M. Elad, "Image denoising via sparse and
%   redundant representations over Learned dictionaries", IEEE Trans.
%   on Image Processing, Vol. 15, no. 12, pp. 3736 - 3745, December 2006.
% 
%   [4] Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli,
%   "Image quality assessment: From error visibility to structural
%   similarity," IEEE Transactios on Image Processing, vol. 13, no. 4,
%   pp. 600 - 612, Apr. 2004.
% 
% 
% Author: Ariel Habshush
% 
% The Cooper Union for the Advancement of Science and Art, 
% Department of Electrical Engineering
% 
% Email: habshu@cooper.edu
% 
% August 2013; Last revision: 15-September-2014
%--------------------------------------------------------------------------
%% Error Checking.
% Check if enough input arguments.
if( nargin < 3)
    error('Not enough input arguments.')
end

% Check if the original image has real valued pixels.
if( any( imag(img(:)) ~=0) )
    error('Image is complexed valued. Please use a real valued image for this function.')
end

% Check if all required parameters exist in the Params struct.
requiredParams = {'PixelRange'; 'numIterCsMRI'; 'lambda'; 'PatchDim';... 
                 'slidingFactor'; 'ompErrorThr'; 'maxNZ';'InitializationMethod';...
                 'adaptDict'; 'noiseSNR'};

tf = isfield(Params,requiredParams);
if( any(tf == 0) )
    missingFields = requiredParams(~tf);
    str1 = 'The following required parameters are missing from the input Params struct:';
    error('%s\n%s\n',str1, missingFields{:} );    
end
clearvars tf str1;

initlMethod = Params.InitializationMethod;
if( strcmpi(initlMethod, 'DCT') || strcmpi(initlMethod,'PCA'))
    if( ~isfield(Params,'numAtoms') )
        error('The method of dictionary initialization requires the number of atoms to be specified.')
    end
elseif( strcmpi(initlMethod, 'External') )
    if( ~isfield(Params,'InitialDict') )
        error('The method of dictionary initialization requires an external dictionary to be provided.')
    elseif( size(Params.InitialDict,1)~= prod(Params.PatchDim))
        error('The number of rows in the external dictionary must equal the number of elements in a patch (i.e. PatchDim(1)*PatchDim(2) )');
    end
else
    strMethods = {'DCT'; 'PCA'; 'External'};
    str1 = 'Invalid dictionary initialization method. Choose from the following methods:';
    error('%s\n%s\n',str1,strMethods{:});
end
clearvars initlMethod strMethods str1;

if( Params.adaptDict == 1)
    if( ~isfield(Params,'numiterateKSVD') || Params.numiterateKSVD <1 )
        error('Specify a valid number of K-SVD iterations');
    end
end
%% Load Parameters from Struct into Variables.
PixelRange = Params.PixelRange;
numIterCsMRI = Params.numIterCsMRI;
lambda = Params.lambda;
PatchDim = Params.PatchDim;
slidingFactor = Params.slidingFactor;
ompErrorThr = Params.ompErrorThr;
maxNZ = Params.maxNZ;
InitializationMethod = Params.InitializationMethod;
adaptDict = Params.adaptDict;
noiseSNR = Params.noiseSNR;
ksvdTrainFactor = Params.ksvdTrainFactor;

% SSIM Paramters
ParamsSSIM.K = Params.SSIM_K;
ParamsSSIM.window = Params.SSIM_window;

if( strcmpi(InitializationMethod, 'DCT') )
    numAtoms = Params.numAtoms;
elseif( strcmpi(InitializationMethod, 'PCA') )
    numAtoms = Params.numAtoms;
elseif( strcmpi(InitializationMethod, 'External') )
    numAtoms = size(Params.InitialDict,2);
end
% Load necessary parameters for K-SVD
if(Params.adaptDict)
ParamsKSVD.displayProgress = 1;
ParamsKSVD.preserveDCAtom = 0; 
ParamsKSVD.errorFlag = 0;
ParamsKSVD.InitializationMethod = 'GivenMatrix';
ParamsKSVD.numIteration = Params.numiterateKSVD;
ParamsKSVD.K = numAtoms;
ParamsKSVD.L = maxNZ;
ParamsKSVD.errorGoal = ompErrorThr;
end
% Notes on above K-SVD parameters: Because all patches have zero mean,
% any patch containing pixels of all the same value will become the
% zero vector. Therefore, there is no requirement to have a constant
% DC atom.
% 
% Miscellaneous Parameters.
% 
%   1) Reduce_DC -- When set to true, all patches will have their mean
%   subtracted out so that they will have zero mean.
Reduce_DC = true; 
%% Prepare Image data for CS-MRI algorithm.
% Convert img to double precision in the range [0.0 - 1.0] for
% processing.
img_orig = img;
img = mat2gray(img_orig,PixelRange);

% Convert noise SNR to a normal distribution standard deviation.
if (noiseSNR ~= Inf)
    sigma_n = sqrt( mean( abs(img(:)).^2 )  ) / ( 10^(noiseSNR/20) );
    img_noisy = addNoise(img, noiseSNR); % Adds Gaussian noise to original image.
else
    sigma_n = 0;
    img_noisy = img;
end

% Compute k-space (Fourier frquency) domain of image.
imgK = fft2(img_noisy);         % Fully sampled k-space of noisy image.
imgK_u = imgK.*kSpaceMask;      % Undersampled k-space data (i.e. the measured data).
img_ZF = ifft2(imgK_u);         % Zero filled k-space Fourier reconstructed image. This image has complex values due to the  IFFT on undersampled k-space.

% Initialize the estimated image to be the magnitude of the  zero
% filled k-space Fourier reconstructed image. Next, break up the
% estimated image into overlapping patches.
img_est = abs(img_ZF); 
[patches, idx] = my_im2col(img_est,PatchDim,slidingFactor);
numPatches = size(patches,2);
numElementsInPatch = size(patches,1);
%% Begin Algorithm
% In order to get better image recovery, the CS-MRI algorithm is ran
% multiple times. Aside from the recovered image and error metrics,
% the following quantities are recored for each iteration:
% 
% 1) numNzInPatch -- Matrix whose j'th column contains the sparsity
% (number of nonzero elements) of the j'th image patch. The i'th row
% corresponds to the i'th algorithm iteration.
% 
% 2) time_total -- Array whose i'th element contains the total time
% (in seconds) it takes to execute the i'th algorithm iteration.
%==========================================================================
% When dealing with dictionary initialization, the DCT and external
% dictionary are fixed so they only have to be declared once (outside
% the iteration loop). On the otherhand, the PCA dictionary is
% dependant on the current estimated image so it has to be created at
% the begining of each iteration.
% 
% Note, the overcomplete DCT dictionary can only have a perfect square
% number of atoms. Therefore, a specified value of numAtoms in the
% Params struct might not be used. Instead, it will be rounded so that
% it will be a perfect square.
%==========================================================================
if( strcmpi(InitializationMethod, 'DCT') )
    Dict_init = overcompleteDCTdict(numAtoms, PatchDim);
    ParamsKSVD.initialDictionary = Dict_init;
    numAtoms = size(Dict_init,2);
elseif( strcmpi(InitializationMethod, 'External') )
    Dict_init = Params.InitialDict;
    ParamsKSVD.initialDictionary = Dict_init;
    numAtoms = size(Dict_init,2);
end

% Preallocate variables.
numNzInPatch = nan(numIterCsMRI, numPatches);
time_total = nan(numIterCsMRI,1);

for iii = 1:numIterCsMRI
    fprintf('csMRI Iteration %d of %d...\n',iii,numIterCsMRI)
    t1 = tic;
    profile on;
    % Subtract mean from patches before sparse coding.
    if(Reduce_DC)
        vecOfMeans = mean(patches);
        patches = patches - repmat(vecOfMeans,numElementsInPatch,1);
    end
    
    numTrainPatch = ksvdTrainFactor*numAtoms;
    if(numPatches >= numTrainPatch)
        idxPatch = randperm( numPatches, numTrainPatch);
    else
        idxPatch = 1:numPatches;
    end
    
    % Initialize PCA dictionary (if specified).
    if( strcmpi(InitializationMethod, 'PCA') )
        Dict_init = pcaDictionary(patches(:,idxPatch), numAtoms );
        ParamsKSVD.initialDictionary = Dict_init;
        if(iii == 1) 
            csMRIout.DictInitial = Dict_init;
        end
    end
    
    % If specified, adapt the dictionary to the estimated image.
    if(adaptDict)
%         [Dict, ksvdOut] = myKSVD(patches(:,idxPatch), ParamsKSVD, ompErrorThr);
        [Dict, ksvdOut] = myKSVD(patches, ParamsKSVD, ompErrorThr);
    else
        Dict = Dict_init;
    end    
    % Perform sparse coding via OMP (or any other sparse coding
    % algorithm) on each of the patches. In order to monitor progress,
    % instead of doing it on all the image patches simultaneously,
    % perform sparse coding on 10,000 patches at a time.
    for jj = 1:10000:numPatches
        jumpSize = min(jj+10000-1, numPatches);
        Coefs = OMPerrn(Dict,patches(:,jj:jumpSize),ompErrorThr,maxNZ);
        numNzInPatch(iii, jj:jumpSize) = sum(Coefs ~= 0);
        
        % Compute patch estimate using the sparse coding coefficients.
        if(Reduce_DC)
            patches(:,jj:jumpSize) = Dict*Coefs + repmat(vecOfMeans(jj:jumpSize),numElementsInPatch,1);
        else
            patches(:,jj:jumpSize) = Dict*Coefs ;
        end
    end
    % Compute the patch-averaged result. For each pixel, add up all of
    % its estimated values from corresponding patches, then divide by
    % the number of patches that pixel was in.
    count = 1;
    numPatchPerPixel= zeros( size(img) );
    img_temp = zeros( size(img) );
    bb = PatchDim(1);
    [rows,cols] = ind2sub( size(img)-bb+1,idx );
    for i  = 1:length(cols)
        col = cols(i); row = rows(i);
        singlePatch = reshape(patches(:,count),[bb,bb]);
        img_temp(row:row+bb-1,col:col+bb-1) = img_temp(row:row+bb-1,col:col+bb-1)+singlePatch;
        numPatchPerPixel(row:row+bb-1,col:col+bb-1) = numPatchPerPixel(row:row+bb-1,col:col+bb-1)+ones(bb);
        count = count+1;
    end;
    img_PatchAvg = img_temp./numPatchPerPixel;
        
    %=========== Apply K-space update formula (Equation (9) in
    %[1])========
    % Compute nu and normalize it by beta, where beta is the number of
    % times each pixel is in a patch (since some pixels are in less
    % patches, we take the maximum).
    %
    % Note: if nu is very large (i.e. noise is very small), then the
    % sampled k-space locations are filled back without taking into
    % account the patch reconstruction update step of our algorithm.
    img_PatchAvgK = fft2(img_PatchAvg);
%     nu = lambda/sigma_n/max(numPatchPerPixel(:)); warning('Nu is
%     inversely proportional to sigma_n and not sigma_n^2');
    nu = lambda/sigma_n^2/max(numPatchPerPixel(:));
    ind = (kSpaceMask==1);
    if(nu>1e20)
        img_PatchAvgK(ind) = imgK(ind);
    else
        img_PatchAvgK(ind) = (1/(1+(nu)))*(img_PatchAvgK(ind) + (nu)*imgK(ind));
    end
    
    % Go from the K-space domain back to the image domain. Then,
    % threshold the estimated image so that its pixel values lie in
    % the range [0.0, 1.0].
    img_est = abs(ifft2(img_PatchAvgK));
    img_est = imThresh( img_est, [0,1]);
    time_total(iii) = toc(t1);
    
    % Compute performance metrics of reconstructed image for each
    % csMRI iteration. The images are converted to original dynamic
    % range when computing the metrics. The function errorMetrics
    % converts it to dynamic range automatically.
    metricsTemp = errorMetrics(img, img_est, PixelRange, ParamsSSIM);    
    metrics.ReconPerIter.MSSIM(iii) = metricsTemp.MSSIM;
    metrics.ReconPerIter.PSNR(iii) = metricsTemp.PSNR;
    metrics.ReconPerIter.SNR(iii) = metricsTemp.SNR;
    metrics.ReconPerIter.MSE(iii) = metricsTemp.MSE;
    metrics.ReconPerIter.RMSE(iii) = metricsTemp.RMSE;
    
    % Break estimated image back into patches for subsequent
    % iteration.
    [patches, idx] = my_im2col(img_est,PatchDim,slidingFactor);
     profile off;
     
end %end csMRI Algorithm

% Store performance metrics of the final reconstructed image and ZF
% Fourier Reconstruction. metrics.Recon = metricsTemp;
metrics.ZF = errorMetrics(img, abs(img_ZF), PixelRange, ParamsSSIM);
metrics.Noisy = errorMetrics(img, img_noisy, PixelRange, ParamsSSIM);
%======================= Not being used
%===================================
% % Threshold the final estimated image so that its pixel values lie
% in the % range [0.0, 1.0] and calculate its error metrics.
% img_est_thr = imThresh(abs(img_est), [0,1]); metrics.ReconThr =
% errorMetrics( abs(img), img_est_thr,maxPixVal );
%==========================================================================

% Return images in a struct.
Images.Original = img;
Images.Recon = img_est;
Images.ZF = abs(img_ZF);
Images.Noisy = img_noisy;


% Return other variables.
csMRIout.DictInitial = Dict_init;
csMRIout.DictFinal = Dict;
csMRIout.NumNzInPatch = numNzInPatch;
csMRIout.AvgNzPerIter = mean(numNzInPatch,2);
csMRIout.sigma_n = sigma_n;
csMRIout.TimePerIteration = time_total;
csMRIout.AvgTimePerIter = mean(time_total);

end %function

