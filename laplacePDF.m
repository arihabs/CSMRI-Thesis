function [ Y ] = laplacePDF( X, Mu, b )
% function [ Y ] = laplacePDF( X, Mu, b )
% 
% Evaluates Laplace distribution at points specified in X.
% 
% Input Parameters:
% 
% X --  Values to evaluate pdf.
% Mu -- Distribution mean.
% b --  Variance = 2b^2.
% 
% Output Parameters:
% Y -- PDF evaluated at X.
%--------------------------------------------------------------------------
Y = ( 1/(2*b) ) * exp( - abs(X - Mu)/b );
end %function

