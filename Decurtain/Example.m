% Removal of curtaining effects by a variational model with directional forward differences
% Author: Jan Henrik Fitschen
% Reference: Jan Henrik Fitschen, Jianwei Ma, Sebastian Schuff, Removal of curtaining effects by a variational model with directional forward differences, Computer Vision and Image Understanding, Volume 155, February 2017, Pages 24-32, ISSN 1077-3142, http://dx.doi.org/10.1016/j.cviu.2016.12.008.
% Email:  jma@hit.edu.cn
% March 22, 2017


load('../MathGeo_testdata/new_artificial.mat');

[ u,s,t ] = curtainDropping( im(:,:,:), 30, .02, .4, .7, 300);

implay(u);