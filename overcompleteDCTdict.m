function [ DCT ] = overcompleteDCTdict( numAtoms, patchDim )
% Generates an overcomplete 2D-DCT dictionary with zero mean atoms.
%
% This function generates the overcomplete dictionary used in "Image
% Denoising Via Sparse and Redundant representations over Learned
% Dictionaries", (appeared in the IEEE Trans. on Image Processing,
% Vol. 15, no. 12, December 2006). This function is based off the code
% provided by the authors.
%--------------------------------------------------------------------------

Pn=ceil(sqrt(numAtoms));
bb = patchDim(1);
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

end

