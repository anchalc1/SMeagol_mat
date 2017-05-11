% xyDet=SM_psf_psfStack_pxInterp(~,xyzEm,opt,reRead)
%
% PSF function using that samples from a given tif stack.
%
% input:
% t         - photon emision times (one per photon, not used here)
% xyzEM     - positions of photon emitters, one xyz triplet per row for
%             each emitted per photon.
% opt       - PSF options, with fields
%   psfTifFile  : path to psf model in the form of a tif stack, should not
%                 be normalized
%   dz          : distance between planes in the psf stack
%   dxy         : in-plane pixel size in psf stack. It is recommended that
%                 these pixels are considerably smaller than the simulated
%                 pixels
%   xyzFocal=[x0 y0 z0] : position of the reference point, e.g., where in
%                         the origin of the model psf model is. x0, y0 are
%                         in pixel units, and z0 must be an integer. 
%                 A good choioce is to take the true position of the
%                 emitter in the focal plane, but other choices (or
%                 localization errors for this point) would only add a
%                 constant offset, so the numbers are not supercritical.
% reRead    - optional flag for re-initializing the persistent variables
%             from the options struct

% output:
% xyDet     - positions of detected photons on the camera chip, randomly
%             generated using the PSF model.
%
%
% M.L. 2017-05-11

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_psfStack_pxInterp, simulated (dummy) PSF in the Smeagol package
% ========================================================================= 
% Copyright (C) 2017 Martin Lindén and Johan Elf
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

%% actual code
function xyDet=SM_psf_psfStack_pxInterp(~,xyzEm,opt,reRead)

if(nargin==0) % return default options struct
    opt=struct;

    opt.dz=50;
    opt.dxy=44;
    opt.xyzFocal=[50 50 31];
    opt.psfTifFile='';

    xyDet=opt;
    
    return;
end
%% construct sampling functions (if not already existing)
persistent pFinite zMin zMax dz dxy x0 y0 z0 NZ z PX PYcX
if(isempty(PX) || (exist('reRead','var') && reRead))
    try
        dxy=opt.dxy;
        dz=opt.dz;
        x0=opt.xyzFocal(1);
        y0=opt.xyzFocal(2);
        z0=opt.xyzFocal(3);
        
        % read the un-normalized psf model
        PSFfile=opt.psfTifFile;
        if(~exist(PSFfile,'file'))
            error([PSFfile ' not found.'])
        end
        pxy=ML_loadStack2(PSFfile);
        
        % compute p(finite) and normalized 2D densities
        NZ=size(pxy,3);
        pFinite=zeros(1,NZ);
        for k=1:NZ
            pFinite(k)=sum(sum(pxy(:,:,k)));
            pxy(:,:,k)=pxy(:,:,k)/pFinite(k);            
        end
        pFinite(k)=pFinite(k)/max(pFinite(k));
        
        z=dz*( (1:NZ)-z0);
        zMin=dz*(1-z0);
        zMax=dz*(NZ-z0);
        
        % pre-compute x-disribution
        PX=cumsum(reshape(sum(pxy,1),size(pxy,2),size(pxy,3)),1);
        
        % pre-compute conditional y-distribution
        PYcX=zeros(size(pxy));
        % normalize
        for k=1:NZ
            for c=1:size(pxy,2)
                zz=sum(pxy(:,c,k));
                if(zz>0)
                    PYcX(:,c,k)=cumsum(pxy(:,c,k))/zz;
                end
            end
        end
        
    catch me
        clear pFinite zMin zMax dz dxy x0 y0 z0 NZ z PX PYcX
        rethrow(me)
    end
end
%% sample
nEm0=size(xyzEm,1); % number of emitted photons
if(nEm0==0)
    xyDet=[];
    return
end
zEm=xyzEm(:,3);

% retain only finite photon positions
withinBounds=zEm>=zMin & zEm<=zMax;
isFinite= withinBounds & rand(nEm0,1)<=interp1(z,pFinite,zEm);
xyzEm=xyzEm(isFinite,:);
nEm=size(xyzEm,1); 
zEm=xyzEm(:,3);

xyDe0=zeros(nEm,2);
xyDet=zeros(nEm,2);
% interpolate z -> psf plane
nz=1+(NZ-1)*(zEm-zMin)/(zMax-zMin);
% sample every point
for m=1:nEm
    % random numbers
    rx=rand;
    ry=rand;
    rdxy=rand(1,2);
    
    % interpolate 2D density
    n0=min(NZ-1,floor(nz(m)));
    dn=nz(m)-n0;
    
    % faster distributions?
    FX=PX(:,n0)*(1-dn)+PX(:,n0+1)*dn;
    xpx=find(FX>=rx,1); %
    
    FY=PYcX(:,xpx,n0)*(1-dn)+PYcX(:,xpx,n0+1)*dn;
    ypx=find(FY>=ry,1);
    xyDet(m,:)=dxy*([xpx-x0 ypx-y0]-0.5+rdxy)+xyzEm(m,1:2);
    
    if(0) % slower, but easier to understand
        % requires pxy to be persistent
        p=pxy(:,:,n0)*(1-dn)+pxy(:,:,n0+1)*dn;
        
        % sample x scatter
        px=sum(p,1);
        px=cumsum(px)/sum(px);
        xpx=find(px>=rx,1); %
        
        % sample conditional y scatter
        py=p(:,xpx);
        py=cumsum(py)/sum(py);
        ypx=find(py>=ry,1);
        
        % convert to length units, and spread evenly over the pixel
        xyDe0(m,:)=dxy*([xpx-x0 ypx-y0]-0.5+rdxy)+xyzEm(m,1:2);
    end
    
end

%num2str(max(abs(xyDe0(:)-xyDet(:))))
dd=sum(withinBounds==0);
if(dd>0)
   warning([ int2str(dd) ' emissions out of z bounds'])
end
