function [ B ] = double2int(A, PixelRange)
% function [ B ] = double2int(A, PixelRange)
% 
% double2int converts a grayscaled image of double precision in the
% range 0.0-1.0 back to its integer representation in the range
% specified by PixelRange.
%--------------------------------------------------------------------------

maxVal = PixelRange(2);
minVal = PixelRange(1);
B = round( A*(maxVal - minVal) + minVal );
end %function

