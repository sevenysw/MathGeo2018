function [A]=OMPerr(D,P,epsilon,options)
%=============================================
% Sparse coding of a group of signals based on a given
% dictionary and specified number of atoms to use.
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  epsilon - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================
options.null = 0;
verbose = getoptions(options,'verbose',0);


[~,L]=size(P);
[n,K]=size(D);
E2 = epsilon^2*n;
% E2 = epsilon*n;
maxNumCoef = n/2;
A = sparse(size(D,2),size(P,2));
if verbose
    tic
end
for k=1:1:L,
    if mod(20*k,20*floor(L/20)) == 0 && verbose
	    time = toc;
        display(['OMPerr Processing ',num2str(floor(100*k/L)),'%, ',num2str(floor(time)),'/',num2str(floor(time/k*L)),' sec'])
    end
    x=P(:,k);
    residual=x;
    indx = [];
    a = [];
    currResNorm2 = sum(residual.^2);
    j = 0;
    while currResNorm2>E2 && j < maxNumCoef,
        j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=D(:,indx(1:j))\x;
        residual=x-D(:,indx(1:j))*a;
        currResNorm2 = sum(residual.^2);
    end;
    if ~isempty(indx)
        A(indx,k)=a;
    end
end;
return;

