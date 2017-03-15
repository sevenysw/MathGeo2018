% This function calculates the frequency parameters of the AGCM model
% 
% Input:
% - AGCMenv, Envelope Parameters from "sepAGCM"
% - data, seismic Data
% - iter, number of Gauss iterations (default: 100)
%
% Output:
% - allPara, all AGCM parameters (envelope and frequency)
% - approx, complete approximation of the data

function [allPara,approx,stepsize] = AGCMfreq_bestStep(EnvPara,data,iter,minStep,k)
multi=2;
% set default values
if nargin<4
    minStep = 10;
end
if nargin<3
    iter=500;
end

% remove envelopes that cause Newton to be very unstable:
% - envelopes with 0 amplitude
rem = abs(EnvPara(1,:))<0.001;
EnvPara = EnvPara(:,~rem);

% calculate number of atoms
NA = size(EnvPara,2);

% define model functions (AGCM envelope + frequency + derivatives)
K = 4;
AGCMenv = @(p,t)(bsxfun(@times,p(1,:).',exp(-bsxfun(@times,p(2,:).',(1-bsxfun(@times,p(3,:).',tanh(K*bsxfun(@minus,t,p(4,:).')))).*bsxfun(@minus,t,p(4,:).').^2))));
AGCMfreq = @(p,t)(cos(bsxfun(@plus,bsxfun(@times,p(5,:).',bsxfun(@minus,t,p(4,:).'))+bsxfun(@times,p(6,:).',bsxfun(@minus,t,p(4,:).').^2),p(7,:).')));
AGCMfreqArg = @(p,t)(bsxfun(@plus,bsxfun(@times,p(5,:).',bsxfun(@minus,t,p(4,:).'))+bsxfun(@times,p(6,:).',bsxfun(@minus,t,p(4,:).').^2),p(7,:).'));
dh = @(p,t)(sin(bsxfun(@plus,bsxfun(@times,p(5,:).',bsxfun(@minus,t,p(4,:).'))+bsxfun(@times,p(6,:).',bsxfun(@minus,t,p(4,:).').^2),p(7,:).')));
%df = @(p,t)(dh(p,t).*bsxfun(@minus,t,p(4,:).'));
%dg = @(p,t)(dh(p,t).*bsxfun(@minus,t,p(4,:).').^2);
dy = @(p,t)(bsxfun(@minus,t,p(4,:).'));

% frequency guess
fdata = fftshift(fft(data));
[~,a] = max(abs(fdata));
freqs = ((-ceil(length(fdata))/2):ceil(length(fdata)/2-1)).';
mF = abs(freqs(a))/length(data)*2*pi;

% phase guess
if k==1; p = log(fdata(a)/abs(fdata(a)))*1i; end
if k==2; p = real(nansum(log(fdata(freqs~=0)./abs(fdata(freqs~=0))).*(1-2*(freqs(freqs~=0)<0))*1i)/length(freqs(freqs~=0))); end
if k==3; p = real(nansum(log(fdata(freqs~=0)./abs(fdata(freqs~=0)))./freqs(freqs~=0)*freqs(a)*1i)/length(freqs(freqs~=0))); end
p = real(p);

% define timesteps and interpolate data to produce more sample points
% (we assume that the data is given such that data(k) is the measured data
% at time k
T = 1:(1/multi):length(data);
data = (interp1(1:length(data),data,T)).';

% take AGCM envelope parameters and add starting values for frequency
% (here very simple: frequency 1, chirp rate 0, phase 0 -> maybe better
% choice?)

allPara = [EnvPara;mF*ones(1,size(EnvPara,2));zeros(1,size(EnvPara,2));p*ones(1,size(EnvPara,2))];
%allPara = [EnvPara;allPara];
%allPara = AGCMstartingParas(EnvPara,data);
%allPara = AGCMstartingParas1(EnvPara,data);

% calculate starting approximation, starting derivative and error
cosMtx = AGCMenv(allPara,T);
sinMtx = cosMtx.*dh(allPara,T);
cosMtx = cosMtx.*AGCMfreq(allPara,T);
approx = sum(cosMtx).';
err(1) = sum((data-approx).^2);
newErr = zeros(minStep,1);

% Gauss iterations
for j=1:iter
    
    % calculate Gauss matrix
    R3 = dy(allPara,T);
    R = sinMtx.*R3;
    R = [R;R.*R3;sinMtx].';
    %R = [AGCMenv(allPara,T).*df(allPara,T);...
    %     AGCMenv(allPara,T).*dg(allPara,T);...
    %     sinMtx].';

    % calculate Gauss step
    dfreq = sparse(R.*(abs(R)>=0.1))\(approx-data);

    % --------------------------------------------------------
    % calculate error for all stepsizes in the following steps
    
    % calculate frequency arguments and cosine+sine matrix for iterations
    FreqArg = 2^(1-minStep)*AGCMfreqArg([allPara(1:4,:);dfreq(1:NA).';dfreq((NA+1):(2*NA)).';dfreq((2*NA+1):end).'],T);
    cosVal = cos(FreqArg);
    sinVal = sin(FreqArg);
    
    % calculate first error
    newErr(end) = sum((data-sum(cosMtx.*cosVal-sinMtx.*sinVal).').^2);
    
    % calculate all other errors
    for k=(minStep-1):(-1):1
    
        % update cosine, sine and frequency values
        sinVal = 2*sinVal.*cosVal;
        cosVal = 2*cosVal.^2-1;
        
        % calculate error
        newErr(k) = sum((data-sum(cosMtx.*cosVal-sinMtx.*sinVal).').^2);
    end
    
    % choose best stepSize
    [~,bestS] = min(newErr);
    stepsize(j) = bestS;
    % check if decreasing
    if newErr(bestS)>err(j)
      %  display(['Stopping in iteration ' num2str(j) ': No decreasing stepsize found!']);
        break;
    end
    
    %update parameters and error
    allPara(5,:) = allPara(5,:) + 2^(1-bestS) * dfreq(1:NA).';
    allPara(6,:) = allPara(6,:) + 2^(1-bestS) * dfreq((NA+1):(2*NA)).';
    allPara(7,:) = allPara(7,:) + 2^(1-bestS) * dfreq((2*NA+1):end).';
    
    % update the amplitudes afterwards
    allPara(1,:) = 1;
    R1 = AGCMenv(allPara,T);
    A = (R1.*AGCMfreq(allPara,T)).';
    allPara(1,:) = sparse(A.*(abs(A)>0.1))\data;

    % update cosine, sine matrix and approximation (maybe better)
    R1 = diag(allPara(1,:))*R1;
    sinMtx = R1.*dh(allPara,T);
    cosMtx = diag(allPara(1,:))*A.';
    approx = sum(cosMtx).';
    % update error
    err(j+1) = sum((data-approx).^2);
    
end


%if j==iter
%    display(['Maximum number of iterations reached!']);
%end

% for return value, reduce number of samples again
approx = approx(1:multi:end);
err = sqrt(err);

end