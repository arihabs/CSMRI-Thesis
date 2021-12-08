function [ img_thr ] = imThresh( img, threshLim )
% function [ img_thr ] = imThresh( img, threshLim )
% 
% Thresholds a real valued image so that pixel values lie in the range
% contained in the 2-element array threshLim.
% 
% Author: Ariel Habshush
% The Cooper Union for the Advancement of Science and Art,
% Department of Electrical Engineering
%
% Email: habshu@cooper.edu
% August 2013; Last revision: 18-July-2014
%--------------------------------------------------------------------------

% Check that the image is real valued.
if( any(imag(img(:))~=0) )
    error('Image threshing not defined for complex valued images.')
end

if( threshLim(1) > threshLim(2) )
    threshLim = [threshLim(2), threshLim(1)];
end

img_thr = img;
img_thr(img < threshLim(1)) = threshLim(1);
img_thr(img > threshLim(2)) = threshLim(2);

end %function

