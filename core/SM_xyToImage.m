% [IM,IMx,IMy]=SM_xyToImage(xy,cam)
% 
% generate image by counting hits of positions xy in the pixel array
% specified by the cam struct: IM(i,j) = counts near pixel center at
% x=IMx(i,j), y=IMy(i,j) , where
% IMx(i,j)=j*cam.pixLength, 1 <= j <= cam.xrange_px
% IMy(i,j)=i*cam.pixLength, 1 <= i <= cam.yrange_px 
%
% Counts further away from cam.pixLength/2 from any pixel center are
% excluded from the image.
% Hence, imagesc(IM) is consistent with 
% plot(xy(:,1)/cam.pixLength,xy(:,2)/cam.pixLength))
%
% cam       - opt.camera from SM options structure. Needs only be passed
%             once, since the function remembers old settings using
%             persistent variables that are updated when this argument is
%             given.
% xy=[x y]  - column vector of positions
%
% IM        - image of pixel counts, for use w imagesc(IM)
% IMx,IMy   - arrays of pixel centers, for use with e.g. surf(IMx,IMy,IM)
%
% M.L. 2014-01-16

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_xyToImage.m, part of the SMeagol package
% ========================================================================= 
% Copyright (C) 2015 Martin Lindén and Johan Elf
% 
% E-mail: bmelinden@gmail.com, johan.elf@gmail.com
% =========================================================================
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.  This program is distributed in the hope that it will
% be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
% the GNU General Public License for more details.
% 
% Additional permission under GNU GPL version 3 section 7
%  
% If you modify this Program, or any covered work, by linking or
% combining it with Matlab or any Matlab toolbox, the licensors of this
% Program grant you additional permission to convey the resulting work.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code
function [IM,IMx,IMy]=SM_xyToImage(xy,cam)

persistent ROI cxy xymin xmax ymax; % update when array changes size
%if( isempty(ROI) || ROI(1)~=cam.xrange_px || ROI(2)~=cam.yrange_px)
if(nargin>1)
    ROI=[cam.xrange_px cam.yrange_px];
    px=cam.pixLength;
    cxy={px*(1:cam.xrange_px),px*(1:cam.yrange_px)};
    xymin=px/2;
    xmax=px*(0.5+ROI(1));
    ymax=px*(0.5+ROI(2));
end

% positions within the camera chip
if(isempty(xy))
    IM=zeros(length(cxy{2}),length(cxy{1}));
else
    %ind=find((xy(:,1)>=xymin).*(xy(:,1)<xmax).*(xy(:,2)>=xymin).*(xy(:,2)<ymax));
    %IM=hist3(xy(ind,:),cxy)'; % bin positions on pixels
    
    IM=hist3(xy(((xy(:,1)>=xymin)&(xy(:,1)<xmax))...
               &((xy(:,2)>=xymin)&(xy(:,2)<ymax)),:),cxy)'; % bin positions on pixels
    % transpose to make consisten coordinate conventions
end

if(nargout>1)
    [IMx,IMy]=meshgrid(cxy{1},cxy{2});
end
