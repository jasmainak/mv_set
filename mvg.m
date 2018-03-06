%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sample from multivariate Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=mvg(N,d,mu,covarm);
% USAGE x=mvg(N,d,mu,covarm)
%
% Generate samples froma multivariate Gaussian distribution with given mean
% and covariance.
%
% INPUT:
%
%   N: Number of samples to generate
%   d: data dimension
%   mu: mean 
%   covarm: covariance matrix 
%
% Example: Generate 100 samples from a 2 dimensional multivariate Gaussian 
% distribution with a mean of [5,6]' and a covariance matrix of [3,0;0,8]:
% x=mvg(100,2,[5,6]',[3,0;0,8])
%

if (nargin~=4),
    error('Not enough input arguments');
end

z=randn(d,N); %generate the samples

[eigVec,eigVal]=eig(covarm);

% apply whitening matrix
%keyboard
x=eigVec*sqrt(eigVal)*z + repmat(mu,[1,N]);
