% Removal of curtaining effects by a variational model with directional forward differences
% Author: Jan Henrik Fitschen
% Reference: Jan Henrik Fitschen, Jianwei Ma, Sebastian Schuff, Removal of curtaining effects by a variational model with directional forward differences, Computer Vision and Image Understanding, Volume 155, February 2017, Pages 24-32, ISSN 1077-3142, http://dx.doi.org/10.1016/j.cviu.2016.12.008.
% Email:  jma@hit.edu.cn
% March 22, 2017

function A = softShrinkage(A,lambda)
%compute soft-shrinkage of A with treshold lambda

if(length(lambda) > 1)
	if(size(A,ndims(A)) ~= length(lambda))
		error('A and lambda must be of same size!');
    end
	
    %adjust lambda to shearlet data structure
	lambda = bsxfun(@times,ones(size(A)),reshape(lambda,[1,1,length(lambda)]));
end

%soft tresholding
aA = abs(A) > lambda;
A = ( A - ( sign(A) .* lambda ) ) .* aA;
