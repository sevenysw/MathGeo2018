%% Asymmetric chirplet transform
% Author: F. Bossmann
% Reference:
% F. Bossmann, J. Ma*, Asymmetric chirplet transform for sparse representation of seismic data, Geophysics, 2015, 80 (6), WD89-WD100.
% F. Bossmann, J. Ma*, Asymmetric chirplet transform ¡ª Part 2: Phase, frequency, and chirp rate, Geophysiscs, 2016, 81(6):V425-V439.
% Email:jma@hit.edu.cn
% April 1st, 2017

%% load and prepare data
display('Performing AGCM transform');
display('loading data...');

% load seismic data matrix
load seismic.mat;
% data = squeeze(DATA).';
% data = data(end:(-1):1,:);
data = data{1}(:,150);
% manipulate all columns to have zero mean and values in [-100,100]
% zero mean is needed for envelope calculation via Hilbert transform
% normalized values allow a general parameter setting
mData = mean(data);
data = bsxfun(@minus,data,mData);
nFactor = 100/max(abs(data(:)));
data = nFactor*data;

% ignore all warnings (some linear systems might be unstable if parameters
% are not chosen wisely)
warning('off','all');

%% envelope reconstruction
fprintf('Calculating envelope parameters...  0%%');
tic;

% define reconstruction parameters
maxRes = 2;                        % stop when infinity norm of residual is less then x% of the data norm (default: 2)
minResChange = norm(data)*0.00025; % stop when the residual changes less than this factor (default: 0.00025*norm(data))
maxIter = 50;                      % number of iterations, also maximal number of AGCM elements (default: 50)

% initialize output parameters
Approx = zeros(size(data));
Parameters = cell(1,size(data,2));

% loop through all columns
for k=1:size(data,2)
    [Approx(:,k),Parameters{k}] = sepAGCM(data(:,k),maxRes,minResChange,maxIter);
    fprintf('%c%c%c%c%3.0f%%',8,8,8,8,k/size(data,2)*100);
end
t = toc;
fprintf('%c%c%c%c',8,8,8,8);
display([' done! (' num2str(t) 's)']);

%% phase reconstruction
fprintf('Calculating phase parameters...  0%%');
tic;

% define reconstruction parameters
phaseChoice = 'all';  % which phase to choose as starting guess (1,2,3 or 'all', default: 'all')
minStep = 10;         % minimal stepsize considered is 2^(-minStep) (default: 10)
maxIter = 500;        % number of maximal iterations (should be dependend on the phaseChoise 1->400, 2->175, 3->325, 'all'->400 or choose 500 slightly better results (performance slowed))

% initialize output parameters

% loop through all columns
if strcmp(phaseChoice,'all')
    for k=1:size(data,2)
        [~,a1] = AGCMfreq_bestStep(Parameters{k},data(:,k),20,minStep,1);
        [~,a2] = AGCMfreq_bestStep(Parameters{k},data(:,k),20,minStep,2);
        [~,a3] = AGCMfreq_bestStep(Parameters{k},data(:,k),20,minStep,3);
        [~,best] = min([norm(a1-data(:,k)),norm(a2-data(:,k)),norm(a3-data(:,k))]);
        [Parameters{k},Approx(:,k),stepsize] = AGCMfreq_bestStep(Parameters{k},data(:,k),maxIter,minStep,best);
        fprintf('%c%c%c%c%3.0f%%',8,8,8,8,k/size(data,2)*100);
    end
else
    for k=1:size(data,2)
        [Parameters{k},Approx(:,k)] = AGCMfreq_bestStep(Parameters{k},data(:,k),maxIter,minStep,phaseChoice);
        fprintf('%c%c%c%c%3d%%',8,8,8,8,k/size(data,2)*100);
    end
end
t = toc;
fprintf('%c%c%c%c',8,8,8,8);
display([' done! (' num2str(t) 's)']);

%% visualization

% count number of atoms used per line
el = zeros(1,size(data,2));
for k=1:size(data,2)
    el(k) = size(Parameters{k},2);
end

display(['MSE: ' num2str(sum((data(:)-Approx(:)).^2)/length(data(:))) ]);
display(['atoms per line: ' num2str(sum(el)/length(el))]);

if size(data,2)>1
    figure;
    subplot(2,2,1); imagesc(data); title('data');
    subplot(2,2,2); imagesc(Approx); title('approximation');
    subplot(2,2,3); imagesc(abs(data-Approx)); title('absolute error');
    subplot(2,2,4); plot(el); title('atoms used');
else
    figure;
    subplot(2,1,1); plot(data); hold on; plot(Approx); title('transform'); legend('data','approx');
    subplot(2,1,2); plot(abs(data-Approx)); title('absolute error');
end

% enable warnings again
warning('on','all');
