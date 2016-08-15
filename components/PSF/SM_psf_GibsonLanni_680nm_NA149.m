% PSF simulation based on an inverse cumulative density function.
% xyDet=SM_psf_GibsonLanni_680nm_NA149(~,xyzEm,opt)
%
% PSF simulation based on an inverse cumulative density function, which is
% given numerically (see below). Up to the maximum Rmax, the output is
% computed as 
% (x_out,y_out) = (x_in,y_in)+ R(u,z_in)*(cos v, sin v),
% where v~U(0,2*pi), and u~U(0,1), i.e., a rotationally symmetric PSF. 
% Additionally, there is a finte z-dependent probability to scatter beyond
% Rmax, in which case the photon is simply discarded (and so it is
% theoretically possible to return xyDet=[]).
%
% This function uses spline interpolation to compute R from a lookup table
% given in a mat-file in
% ../PSFdata/numericalPSF_GibsonLanni_584nm_NA14.mat, which in turn was
% created using results from the PSFgenerator software [1-3], with standard
% settings for the Gibson-Lanni PSF model[4], 680 nm, and NA=1.49, and 5 nm
% voxel size in the PSF simulations.
% 
% Photons outside the parameterized z-range (-+ 6350) are treated using the
% nearest z value, and a warning is issued. 
%
% input:
% t         - photon emision times (one per photon, not used here)
% xyzEM     - positions of photon emitters, one xyz triplet per row for
%             each emitted per photon.
% opt       - PSF options (not used here)
%
% output:
% xyDet     - 1) if no input arguments are given, an empty struct (default
%                options) are returned. 
%             2) if all input argument are given, xyDet are positions of
%             detected photons on the camera chip, randomly generated using
%             the PSF model.
%
% Parameterize you own: numericalPSF_GibsonLanni_680nm_NA149.mat has fields
% u,z,R_uz,PfiniteR, which parameterizes the function R(u,z_in) and the
% z-dependent probability to discard a photon. 
% Sizes:    R_uz~nu*nz, u~1*nu, z~1*nz, PfiniteR~1*nz, so that 
%           [z,u]=meshgrid(D.z,D.u);
%           surf(D.R_uz,z,u,'edgecolor','none')
%           looks nice.
% -------------------------------------------------------------------------
% NOTE : This way to parameterize a PSF can lead to artificially high
% intensities near the center of the spot (R=0). This should be OK in
% practise as long as the voxel size used in parameterization (5 nm) is
% significantly smaller than the pixels used by SMeagol, but no
% quantitative guarantees are given...
% -------------------------------------------------------------------------
% refs: 
% 1. http://bigwww.epfl.ch/algorithms/psfgenerator/
% 2. Hagai Kirshner, François Aguet, Daniel Sage, Michael Unser, 3-D PSF
% Fitting for Fluorescence Microscopy: Implementation and Localization
% Application, Journal of Microscopy, vol. 249, no. 1, pp. 13-25, January
% 2013.   
% 3. Alessandra Griffa, Nathalie Garin and Daniel Sage, Comparison of
% Deconvolution Software in 3D Microscopy. A User Point of View, Part I and
% Part II, G.I.T. Imaging & Microscopy, vol 1, pp. 43-45, 2010
% 4. Gibson, S. F. & Lanni, F. Experimental test of an analytical model of
% aberration in an oil-immersion objective lens used in three-dimensional
% light microscopy. J Opt Soc Am A 9, 154–166 (1992).
% ------------------------------------------------------------------------
% M.L. 2015-01-19

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_GibsonLanni_680nm_NA149, PSF-simulator in the SMeagol package
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

%% actual code
function xyDet=SM_psf_GibsonLanni_680nm_NA149(~,xyzEm,~)
if(nargin==0) % return default options struct
    opt=struct();
    xyDet=opt;
    return;
end
%% construct sampling functions (if not already existing)
persistent FuMax FR_uz zMin zMax;
if(isempty(FuMax))
    try
        PSFfile=fullfile(mfilename('fullpath'),'..','..','PSFdata','numericalPSF_GibsonLanni_680nm_NA149.mat');
        if(~exist(PSFfile,'file'))
            error('SM_psf_GibsonLanni_680nm_NA149: PSF data file numericalPSF_GibsonLanni_680nm_NA149.mat not found.')
        end
        data=load(PSFfile,'u','z','R_uz','PfiniteR');
        z=data.z;
        u=data.u;
        R_uz=data.R_uz;
        PfiniteR=data.PfiniteR;
        
        % construct interpolants
        FR_uz=griddedInterpolant({u,z},R_uz,'spline','nearest');         % inverse radial cdf
        FuMax=griddedInterpolant(z,PfiniteR,'spline','nearest');
        zMin=min(z);
        zMax=max(z);
    catch me
        clear FR_uz FuMax zMin zMax
        rethrow(me)
    end
end
%% sample radially symmetric scattering events
xyDet=[]; % default answer if no photons scatter to finite radius
if(size(xyzEm,1)>0)
    
    indW=find((xyzEm(:,3)<zMin)+(xyzEm(:,3)>zMax),1);
    if(~isempty(indW))
        warning('SM_psf_GibsonLanni_680nm_NA149: input z out of range')
    end
    
    u=rand(size(xyzEm,1),1); % one U(0,1) random number for each photon
    uMax=FuMax(xyzEm(:,3));  % probability to get a finite scattering radius
    
    indF=find(u<=uMax);      % photons to get a finite position
    if(~isempty(indF))
        
        Rout=zeros(length(indF),1); % only keep photons with finite scattering radius
        u=rand(size(Rout));         % new random numbers to determine the scattering radius
        Rout=FR_uz(u,xyzEm(indF,3));% scattering radius
        
        theta=2*pi*rand(size(u)); % scattering angle: uniformly distributed
        xyDet=xyzEm(indF,1:2)+[Rout.*cos(theta) Rout.*sin(theta)];
    end
end
end
