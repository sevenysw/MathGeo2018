addpath(genpath('./'))

if ~exist('readsegy','file')
    disp('download Seismiclab');
    untar('http://seismic-lab.physics.ualberta.ca/SeismicLab.tar.gz','.');
    disp('include Seismiclab');
    addpath(genpath('SeismicLab'));
end