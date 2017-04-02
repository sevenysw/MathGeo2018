function  Dc  = clip(D,cmin,cmax);
%CLIP: A program to clip seismic data.
%
%  Dc  = clip(D, cmin, cmax);
%
%  IN   D:          data to be clipped
%       cmin:       lower clip in % (90% means clip at the 90% negative amplitude)
%       cmax:       lower clip in % (90% means clip at the 90% positive amplitude)
%
%  OUT  Dc:         data after being clipped
%
%  Example: 
%
%  d = sin(2*pi*.02*[1:1:200]); plot(clip(d,90,90));
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%


 if (nargin < 1 | nargin > 3)
  error('Wrong number of input parameters in CLIP.')
 end

%   [nt,nx] = size(D);

%   d = reshape(D,nx*nt,1);
  d = D(:);

  dmin = abs(min(d));
  dmax = abs(max(d));
        
  I = find(d>0 & d > cmax*dmax/100.);
  J = find(d<0 & abs(d) > cmin*dmin/100);

  d(I) = cmax*dmax/100.;
  d(J) = -cmin*dmin/100.;

%   Dc = reshape(d,nt,nx);
  Dc = reshape(d,size(D));
  return;
