%% Asymmetric chirplet transform
% Author: F. Bossmann
% Reference:
% F. Bossmann, J. Ma*, Asymmetric chirplet transform for sparse representation of seismic data, Geophysics, 2015, 80 (6), WD89-WD100.
% F. Bossmann, J. Ma*, Asymmetric chirplet transform ¡ª Part 2: Phase, frequency, and chirp rate, Geophysiscs, 2016, 81(6):V425-V439.
% Email:jma@hit.edu.cn

% Function calculates the AGCM envelope approximation of Data
% Input:
%   Data = input data
%   maxRes = maximal resdual value as stopping criterion
%   minResChange = minimal resdual norm change as stopping criterion
%   maxIter = maximal number of iterations as stopping criterion
%   eps = epsilon for choosing the interval Ifist:Ilast
%
% Output:
%   approx = approximated envelope
%   para = used AGCM parameters
%   res = residual
%
% Note: The input data "Data" should be a vector with mean 0!
% April 1st, 2017

function [approx,para,res] = sepAGCM(Data,maxRes,minResChange,maxIter,eps)

% set input parameters standard values
if nargin<5 || isempty(eps)
    eps = 0.2;
end
if nargin<4 || isempty(maxIter)
    maxIter = 100;
end
if nargin<3 || isempty(minResChange)
    minResChange = norm(Data)*0.00025;
end
if nargin<2 || isempty(maxRes)
    maxRes = 2;
end

% calculate envelope
Data = abs(hilbert(Data));


% define AGCM envelope
K = 4;
AGCMenv = @(p,t)(p(1)*exp(-p(2)*(1-p(3)*tanh(K*(t-p(4)))).*(t-p(4)).^2));

% set starting values
approx = zeros(size(Data));
res = Data;
lastResNorm = inf;
para = zeros(4,maxIter);
A = [];

T = (1:length(approx)).';

% start MP loop
while   (max(abs(res))>maxRes) && ...
        ((lastResNorm-norm(res))>minResChange) && ...
        maxIter>0
    
    % starting guess
    p = zeros(4,1);
    [p(1),p(4)] = max(res); % starting guess amplitude and time shift
    
    % find first index that is above threshold, i.e. starting point of
    % interval I
    firstI = find(res(1:p(4))<(abs(p(1))*eps),1,'last')+1;
    if isempty(firstI)
        firstI = 1;
    end
    if any(diff(res(firstI:p(4))))<-eps
        firstI = firstI+find(diff(res(firstI:p(4)))<-eps,1,'last');
    end
        
    % find last index
    lastI = find(res(p(4):end)<(abs(p(1))*eps),1,'first')+p(4)-2;
    if isempty(lastI)
        lastI = length(res);
    end
    if any(diff(res(p(4):lastI))>eps)
        lastI = p(4)-1+find(diff(res(p(4):lastI))>eps,1,'first');
    end
    I = (firstI:lastI).';
    
    % left decay approximate
    ts = (firstI:(p(4)-1)).';
    if isempty(ts)
        lD = 0;
    else
        lD = median(log(abs(p(1))./res(ts))./((ts-p(4)).^2));
    end
    
    % right decay approximate
    ts = ((p(4)+1):lastI).';
    if isempty(ts)
        rD = 0;
    else
        rD = median(log(abs(p(1))./res(ts))./((ts-p(4)).^2));
    end

    p(2) = (lD+rD)/2;
    p(3) = (lD-rD)/(lD+rD);
    if isnan(p(3))
        p(3) = 0;
    end
    
    % update
    para(1:4,maxIter) = p;
    lastResNorm = norm(res);
    
    % optimize parameters
    A = [AGCMenv([1;para(2:4,maxIter)],T),A];
    relI = sum(abs(A),2)>0.2;
    para(1,maxIter:end) = lsqnonneg(A(relI,:),Data(relI));

    approx = A*para(1,maxIter:end).';
    res = Data-approx;
    maxIter = maxIter-1;
end

% prepare output
para = para(:,(maxIter+1):end);

end
