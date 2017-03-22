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
