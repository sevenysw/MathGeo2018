function Iout=InpaintingInterp2(Ip,Fmask,interp)
%Pre-interpolation with matlab function griddata.

%Ip:    subsampled data
%Fmask: sampling mask
%interp:interpolation method

%Iout: interpolated data


if (~exist('interp','var'))
    interp='nearest';
end

Ip=double(Ip);
[N,M]=size(Ip);
[X,Y]=meshgrid([1:M],[1:N]);
x=X(Fmask==1);
y=Y(Fmask==1);
z=Ip(Fmask==1);

xi=X(Fmask==0);
yi=Y(Fmask==0);


Iout=Ip;
interp=lower(interp);
zi = griddata(x,y,z,xi,yi,interp);
zi(isnan(zi))=0.;
Iout(Fmask==0)=zi;
