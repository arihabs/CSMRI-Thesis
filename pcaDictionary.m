function Dict = pcaDictionary(patches, numAtoms )
% function Dict = pcaDictionary(patches, numAtoms )
% 
% pcaDictionary.m Generates a sparsifying dictionary using principle
% component analysis (PCA). This code is largely based on the code
% used in the DLMRI software package used in [1].
% 
% INPUT PARAMETERS:
% 
%       1) patches -- matrix whose i'th column contains the i'th image
%       patch (arranged in lexicographical order) used to train the
%       PCA dictionary.
% 
%       2) numAtoms -- desired number of atoms (i.e. rows) in the
%       resulting PCA dictionary. If the numAtoms is greater than the
%       number of training patches, columns 1 through numAtoms of PCA
%       dictionary contain the right singular values while the
%       remaining atoms are filled in with a random subset of training
%       patches. If the numAtoms is less than the number of elements
%       in a single training patch, columns 1 through numAtoms of PCA
%       dictionary contain the right singular values corresponding to
%       the "numAtoms" greatest singular values.
% 
% Code Description: PCA is performed on the image patches to form a
% sparsiying dictionary. By taking the singular value decomposition
% (SVD) of the collection of image patches, the right singular vectors
% form an orthogonal basis such that the patch data has the greatest
% variance in the direction of the first vector, the second greatest
% variance the direction of the second vector, and so on. This allows
% for the possibilty of representing the data as a K-sparse vector if
% the data is transformed using the first K singular vectors.
% 
% In order to speed up the svd computation, we can multiply
% patches*patches' which gives us a smaller matrix (assuming that the
% number of elements in a patch is smaller than the number of patches
% we are training on) to perform svd on. The left singular vectors of
% this smaller matrix is equivalent to the right singuler vector of
% the original patches.
% 
% NOTE: In general, when performing performing PCA on any data matrix
% X, X takes the form such that its rows represent different
% experiment measurements and its columns represent different
% features. This function assumes the data matrix is of the transpose
% form such that the columns represent different measurements (i.e.
% different patches) and the rows represent different features (i.e.
% different pixel locations).
% 
% 
% REFERENCES:
% 
%   [1] S.Ravishankarand Y. Bresler, "MR image reconstruction from
%   highly undersampled k-space data by dictionary learning," IEEE
%   Trans. Med. Imag., vol. 30, no. 5, pp. 1028 - 1041, 2011.
% 
% 
% Author: Ariel Habshush
% 
% The Cooper Union for the Advancement of Science and Art,
% Department of Electrical Engineering
%
% Email: habshu@cooper.edu
% 
% August 2013; Last revision: 20-July-2014
%--------------------------------------------------------------------------

numElementsInPatch = size(patches,1);
numPatches = size(patches,2);
Dict = nan(numElementsInPatch, numAtoms);
[U,~,~] = svd(patches*patches');
if (numAtoms < numElementsInPatch)
    Dict(:,1:numAtoms) = U(1:numAtoms);
else
    Dict(:,1:numElementsInPatch) = U;
end
numRemainAtoms = numAtoms-numElementsInPatch;
if( numRemainAtoms > 0 )
    idxPatch = randperm(numPatches,numRemainAtoms);
    Dict(:,numElementsInPatch+1:end) = patches(:,idxPatch);
end
        
        
        
end %function

