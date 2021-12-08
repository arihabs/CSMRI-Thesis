% Script used to learn an image prior via the K-SVD algorithm.
%--------------------------------------------------------------------------

%% Add K-SVD Functions to Path
addpath(genpath('./KSVD/'));
presentDate =  datestr(now,'yyyymmdd');
%% Randomly choose patients and slices for training and testing.
randselect = false;
load('subject_and_slice_ranges.mat');
numSubjects = length(subjectRange);
numSlices = length(sliceRange);
if(randselect)
    numTrainSubjects = 15;
    trainSubjects = subjectRange(randperm(numSubjects,numTrainSubjects));
    valSubjects = setdiff(subjectRange,trainSubjects);
else
    load('training_validation_set.mat');
    numTrainSubjects = length(trainSubjects);
end
%% Initialize K-SVD Parameters.
pctSparse = 0.125;
redundacyFactor = 2;
ParamsKSVD.errorGoal = 0;
ParamsKSVD.numIteration = 180;
ParamsKSVD.PatchDim = [8,8];
numElementsPerPatch = prod(ParamsKSVD.PatchDim);
numAtoms = redundacyFactor * numElementsPerPatch;
Dict_init = overcompleteDCTdict(numAtoms, ParamsKSVD.PatchDim );
numAtoms = size(Dict_init,2);
ParamsKSVD.initialDictionary = Dict_init;

ParamsKSVD.K = numAtoms;
ParamsKSVD.L = round( pctSparse*numElementsPerPatch );
ParamsKSVD.displayProgress = 1;
ParamsKSVD.preserveDCAtom = 0; 
ParamsKSVD.errorFlag = 0;
ParamsKSVD.InitializationMethod = 'GivenMatrix';
ParamsKSVD.initialDictionary = Dict_init;
%% Gather patches from training subjects for training
varSort = true; % set true if patches with highest variance should be chosen as training signals.
numTrainingSignals = 50e3; % Number of training signals (patches) to use in KSVD.
numSlicePerSubject = 2;
slidingFactor = 2;
numSignalsPerSubject = numTrainingSignals/numTrainSubjects;
numSignalsPerSlice = round(numSignalsPerSubject/numSlicePerSubject);

filelocation = '..\Datasets\TwentyNormalBrains\Entire_Brain'; % chanage
PixelRange = [0,2^12-1];
% Keep track of slices chosen.
sliceNum = zeros(numTrainSubjects,numSlicePerSubject);

trainingPatches = [];
for itrSubjct = 1:numTrainSubjects % Loop through subjects
    currentSubject = trainSubjects(itrSubjct);
    sliceIdx = randperm( numSlices, numSlicePerSubject ); % choose random slices
    sliceNum(itrSubjct,:) = sliceRange(sliceIdx);
    
    filename = sprintf('subject%02d_t1w_p4.mat',currentSubject);
    load( fullfile(filelocation,filename),'imgVol');    
    for itrSlice = sliceIdx % Loop through slices
        currentSlice = imgVol(:,:,itrSlice);
        [patches, idx] = my_im2col( mat2gray(currentSlice,PixelRange), ParamsKSVD.PatchDim, slidingFactor);
        numPatches = size(patches,2);
        if(varSort)
            [~, varIdx] = sort(var(patches),'descend');
            patchIdx = varIdx( 1:min(numSignalsPerSlice,numPatches) ); %choose patches with largest variance in order to exclude black background from training.       
        else
            patchIdx = randperm(numPatches, min(numSignalsPerSlice, numPatches) );
        end
        trainingPatches = [trainingPatches, patches(:,patchIdx)];
    end
end
% figure; montage3Darray( reshape(trainingPatches,[ParamsKSVD.PatchDim,size(trainingPatches,2)] ) )

% Remove DC !
vecOfMeans = mean(trainingPatches);
trainingPatchesNoDC = trainingPatches - repmat(vecOfMeans,numElementsPerPatch,1);
%% Perform K-SVD on Training Patches.
t1 = tic;
[ Dict, ksvdOut ] = myKSVD(trainingPatchesNoDC, ParamsKSVD, ParamsKSVD.errorGoal);
t_total = toc(t1);
%% Save Results.
DictParams.numTrainSignals = numTrainingSignals; 
DictParams.slidingFactor = slidingFactor;
% DictParams.pctKspaceSampled = pctSampled;
DictParams.numKSVDIter = ParamsKSVD.numIteration;
DictParams.sparsity = ParamsKSVD.L;
DictParams.patchDim = ParamsKSVD.PatchDim;
DictParams.ompErrorGoal = ParamsKSVD.errorGoal;
DictParams.initialDict  = '2D-DCT';

DictResults.ksvdOut = ksvdOut;
DictResults.avgTimePerIter = mean(ksvdOut.TimePerIteration);
DictResults.totalTime = datestr(t_total/24/3600, 'HH:MM:SS');
DictResults.slicesUsed = [trainSubjects',sliceNum];
DictResults.numAtoms = size(Dict,2);

GlobalDict.dict = Dict;
GlobalDict.params = DictParams;
GlobalDict.results = DictResults;


