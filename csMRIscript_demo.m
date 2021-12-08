% Demo of how to use csMRI.m function.
%
% This script performs CS-MRI on 3 brain scan images and recovers the
% image using 2D-DCT, DLMRI (here refered to as "PCA"), and Prior
% dictionary (also refered to as "Global" in this code). Averages of
% different metrics are then computed across the subjects.
%
% *NOTE 1: Users must change the following variables to conform with
% their setup:
%
%   dictPath -- relative or absolute path where prior dictionary is
%   located.
% 
%   dictFile -- name of prior dictionary file (including .mat
%   extension).
% 
%   filelocation -- location containing brain scans (in .mat format).
% 
%   filename -- brain scan file name.
% 
%   testSet -- the subject number and slice number must be valid for
%   the dataset being used. If the BRAINWEB dataset is used, then the
%   values here are valid.
%
% *NOTE 2: Make sure all functions called by csMRI.m are added to the
% path.
%
% All results are stored in the multidimensional struct (each
% dimension corresponding to a test subject) called "Subject". The
% struct contains fields for each type of recovery method (DCT, PCA,
% Prior). Each of those fields contain subfields with the results of
% csMRI.m. See csMRI.m for information about those fields.
%--------------------------------------------------------------------------

%% Add to path all required m-files.
addpath(genpath('./KSVD/'));
%addpath('../MATLABpackages/ssim/');
%% Specify test subject and slice number to perform CS-MRI.
% First column of testSet contains the subject number and the second
% column contains the slice index.
testSet = [20,12;
           38,34;
           43,10];
numValSbjcts = size(testSet,1);
%% Specify amount of compression and generate sampling mask.
imgDim = [256,256];  % Image dimension.
pctSampled = 0.40;   % Percentage of k-space to sample.
pctCenter = 0.5;     % Percentage points to be taken from the center of k-space.
pdfType = 'Normal';  % Sample from a Gaussian distribution.
pdfParam.Mu = floor(imgDim(1)/2);
pdfParam.Sigma = 70;
SampleMask = genSamplingMask(imgDim,pctSampled,pctCenter,pdfType, pdfParam);
%% Set CS-MRI Parameters.
Params.numAtoms = 144;
Params.noiseSNR = 15;
Params.maxNZ = 8;
Params.PatchDim = [8,8];
numElementsPerPatch = prod(Params.PatchDim);
Params.ompErrorThr = 0.025;
Params.ksvdTrainFactor = 200;
Params.numIterCsMRI = 15;
Params.numiterateKSVD = 10;
Params.PixelRange = [0,2^12-1];
Params.lambda = 0.1;
Params.C = 0;
Params.slidingFactor = 1;
Params.SSIM_K = [0.05,0.05];
Params.SSIM_window = fspecial('gaussian', 11, 1.5);
%% Load Prior Dictionary
% Specify the prior dictionary file path and the file name.
dictPath  = sprintf( './GlobalDicts./Entire_Images./%dX%d',...
                         Params.PatchDim(1),Params.PatchDim(2) );

dictFile = sprintf('GlobalDictKSVD_%dX%dpatch_%dAtoms.mat',...
                    Params.PatchDim(1),Params.PatchDim(2),Params.numAtoms); 
load( fullfile(dictPath,dictFile));
%% Iterate over each Subject.
field1 = 'SubjectNumber';
field2 = 'SliceNumber';
Subject = struct( field1, num2cell(testSet(:,1)), field2, num2cell(testSet(:,2)) );

for itrSbjct = 1:numValSbjcts
    %% Load a single image slice from the volume data.
    currentSbjct = testSet(itrSbjct,1);
    currentSliceIdx = testSet(itrSbjct,2);    
    filelocation = '.\Datasets\TwentyNormalBrains\Entire_Brain';
    filename = sprintf('subject%02d_t1w_p4.mat',currentSbjct);
    load( fullfile(filelocation,filename),'imgVol');
    img_orig = imgVol(:,:,currentSliceIdx);
    clearvars imgVol;
        
    %% Run CSMRI for different Dictionaries.
    fprintf('Processing Subject %d of %d... \n',itrSbjct, numValSbjcts);
    disp('Processing DCT...');
    
    % 2D-DCT.
    ParamsDCT = Params;
    ParamsDCT.InitializationMethod = 'DCT'; % Choose from 'DCT', 'PCA', 'External'
    ParamsDCT.adaptDict = 0;
    ParamsDCT.InitialDict = [];
    [Images, Metrics, CSMRIout] = csMRI(img_orig,SampleMask.Mask,ParamsDCT);
    Subject(itrSbjct).DCT.Images = Images;
    Subject(itrSbjct).DCT.Metrics = Metrics;
    Subject(itrSbjct).DCT.CSMRIout = CSMRIout;

    disp('DCT Complete!');

    % PCA (DLMRI Method)
    disp('Processing PCA...');

    ParamsPCA = Params;
    ParamsPCA.InitializationMethod = 'PCA'; % Choose from 'DCT', 'PCA', 'External'
    ParamsPCA.adaptDict = 1;
    ParamsPCA.InitialDict = [];
    [Images, Metrics, CSMRIout] = csMRI(img_orig,SampleMask.Mask,ParamsPCA);

    Subject(itrSbjct).PCA.Images = Images;
    Subject(itrSbjct).PCA.Metrics = Metrics;
    Subject(itrSbjct).PCA.CSMRIout = CSMRIout;

    disp('PCA Complete!');

    % Prior Dictionary
    disp('Processing Prior Dictionary...');
    ParamsGlobal = Params;
    ParamsGlobal.InitializationMethod = 'External'; % Choose from 'DCT', 'PCA', 'External'
    ParamsGlobal.adaptDict = 0;
    ParamsGlobal.InitialDict = GlobalDict.dict;
    [Images, Metrics, CSMRIout] = csMRI(img_orig,SampleMask.Mask,ParamsGlobal);

    Subject(itrSbjct).Global.Images = Images;
    Subject(itrSbjct).Global.Metrics = Metrics;
    Subject(itrSbjct).Global.CSMRIout = CSMRIout;

    disp('Global Prior Complete!');
end
%% Calculage Averages over Subjects.
Recon_PSNR = nan(numValSbjcts, Params.numIterCsMRI,3);
Recon_MSSIM = nan(numValSbjcts, Params.numIterCsMRI,3);
Noisy_PSNR = nan(numValSbjcts,3);
Noisy_MSSIM = nan(numValSbjcts,3);
ZF_PSNR = nan(numValSbjcts,3);
ZF_MSSIM = nan(numValSbjcts,3);
Recon_AvgTimePerIter = nan(numValSbjcts,3);

dctDim = 1;
pcaDim = 2;
priorDim = 3;

for itrSbjct = 1:numValSbjcts
    Recon_PSNR(itrSbjct,:,dctDim) = Subject(itrSbjct).DCT.Metrics.ReconPerIter.PSNR;
    Recon_PSNR(itrSbjct,:,pcaDim) = Subject(itrSbjct).PCA.Metrics.ReconPerIter.PSNR;
    Recon_PSNR(itrSbjct,:,priorDim) = Subject(itrSbjct).Global.Metrics.ReconPerIter.PSNR;
    Recon_MSSIM(itrSbjct,:,dctDim) = Subject(itrSbjct).DCT.Metrics.ReconPerIter.MSSIM;
    Recon_MSSIM(itrSbjct,:,pcaDim) = Subject(itrSbjct).PCA.Metrics.ReconPerIter.MSSIM;
    Recon_MSSIM(itrSbjct,:,priorDim) = Subject(itrSbjct).Global.Metrics.ReconPerIter.MSSIM;

    Noisy_PSNR(itrSbjct,dctDim) = Subject(itrSbjct).DCT.Metrics.Noisy.PSNR;
    Noisy_PSNR(itrSbjct,pcaDim) = Subject(itrSbjct).PCA.Metrics.Noisy.PSNR;
    Noisy_PSNR(itrSbjct,priorDim) = Subject(itrSbjct).Global.Metrics.Noisy.PSNR;
    Noisy_MSSIM(itrSbjct,dctDim) = Subject(itrSbjct).DCT.Metrics.Noisy.MSSIM;
    Noisy_MSSIM(itrSbjct,pcaDim) = Subject(itrSbjct).PCA.Metrics.Noisy.MSSIM;
    Noisy_MSSIM(itrSbjct,priorDim) = Subject(itrSbjct).Global.Metrics.Noisy.MSSIM;

    ZF_PSNR(itrSbjct,dctDim) = Subject(itrSbjct).DCT.Metrics.ZF.PSNR;
    ZF_PSNR(itrSbjct,pcaDim) = Subject(itrSbjct).PCA.Metrics.ZF.PSNR;
    ZF_PSNR(itrSbjct,priorDim) = Subject(itrSbjct).Global.Metrics.ZF.PSNR;
    ZF_MSSIM(itrSbjct,dctDim) = Subject(itrSbjct).DCT.Metrics.ZF.MSSIM;
    ZF_MSSIM(itrSbjct,pcaDim) = Subject(itrSbjct).PCA.Metrics.ZF.MSSIM;
    ZF_MSSIM(itrSbjct,priorDim) = Subject(itrSbjct).Global.Metrics.ZF.MSSIM;

    Recon_AvgTimePerIter(itrSbjct,dctDim) = Subject(itrSbjct).DCT.CSMRIout.AvgTimePerIter;
    Recon_AvgTimePerIter(itrSbjct,pcaDim) = Subject(itrSbjct).PCA.CSMRIout.AvgTimePerIter;
    Recon_AvgTimePerIter(itrSbjct,priorDim) = Subject(itrSbjct).Global.CSMRIout.AvgTimePerIter;

end

% Calculate Mean+stdPSNR and mean+std MSSIM.
Recon_PSNRavg = mean( squeeze(Recon_PSNR(:,end,:)) );
Recon_PSNRstd = std( squeeze(Recon_PSNR(:,end,:)) );
Recon_MSSIMavg = mean( squeeze(Recon_MSSIM(:,end,:)) );
Recon_MSSIMstd = std( squeeze(Recon_MSSIM(:,end,:)) );

Noisy_PSNRavg = mean( Noisy_PSNR(:) );
Noisy_PSNRstd = std( Noisy_PSNR(:) );
Noisy_MSSIMavg = mean( Noisy_MSSIM(:) );
Noisy_MSSIMstd = std( Noisy_MSSIM(:) );

ZF_PSNRavg = mean( ZF_PSNR(:) );
ZF_PSNRstd = std( ZF_PSNR(:) );
ZF_MSSIMavg = mean( ZF_MSSIM(:) );
ZF_MSSIMstd = std( ZF_MSSIM(:) );

Recon_AvgTime = mean(Recon_AvgTimePerIter,1);
%% Sample Plot
figure; 
imshow(Subject(1).Global.Images.Original);
title('Original Image');

figure; 
imshow(Subject(1).Global.Images.ZF);
title('Zero-Filled k-space Reconstruction');

figure; 
imshow(Subject(1).Global.Images.Recon);
title('Prior Dictionary Reconstruction');