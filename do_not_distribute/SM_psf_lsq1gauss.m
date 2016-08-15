% function [xyDet,sigma]=SM_psf_lsq1gauss(t,xyzEm,opt)
%
% Explicit PSF using a least-squares fit of a single Gaussian to images of
% 40 nm beads, taken at 561 nm with Andor Ixon Ultra (Elf lab experiment
% EXP-14-BH0518).
%
% The model is the Gaussian standard deviation as a function of z, fitted
% in the approximate range ~[-1250 , 375]  nm. sig ~< 600 nm in this range
% (for 561 nm wavelength), and the minimum sigma is for z=0. z values
% outside this range are changed to the closest value inside the range, and
% a warning is issued. For wavelengths other than 561 nm, the std deviation
% is rescaled as sig(z,L)= L/561nm*sigma(z,561nm), but the z-dependence is
% not rescaled.
%
% NOTE: the PSF-width depends strongly and asymmetrically on z, and so it
% is provbably good to have more sample volume below the focal plane (z=0)
% than above it. A plot of sigma vs z should give some guidance in how to
% place the focal plane.
%
% output:
% xyDet     - positions of detected photons on the camera chip, randomly
%             generated using the PSF model.
% sigma     - PSF model std for every emitted photon, which depends on the
%             simulated wave-length and z-positions in xyzEm(:,3). This
%             output is useful for checking the output. 
% input:
% t         - photon emision times (one per photon, not used here)
% xyzEM     - positions of photon emitters, one xyz triplet per row for
%             each emitted per photon.
% opt       - PSF options
% opt.lambda- light wavelength in simulation (L in the above equation)
% 
%
% The options object opt here corresponds to the opt.psf structure defined
% in the Palantir runinput file. Calling SM_psf_lsq1gauss() with no
% parameters returns an options object with default parameters.
%
% M.L. 2015-01-15

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_lsq1gauss, simulated PSF in the SMeagol package
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
function [xyDet,sig]=SM_psf_lsq1gauss(varargin)

if(nargin==3) % prepare to do stuff
    xyzEm=varargin{2};
    opt=varargin{3};    
elseif(nargin==0) % return a struct with default parameter values
    opt.lambda=561; % nm    
    xyDet=opt;
    sig=[];
    return
else % something wrong
   error('SM_psf_lsq1gauss requires three (or zero) input parameters')
end


% model parameters
Psigma=[3.085758651e-18  9.562807673e-15  1.161878277e-11   6.919807998e-09 1.963969739e-06  0.0006148163666 -6.113601909e-08      137.6250611];

lambda0=561; % wave-length at which data was taken
dzMin=-1250; % interval where the polynomial fit is valid
dzMax=375;

lambda1=opt.lambda; % laser wave length for this simulation

% emitted photons
Nph=size(xyzEm,1); % number of emitted photons

% out-of-focus heights
dz=xyzEm(:,3);

% check for out-of-bounds
if(~isempty(find(dz<dzMin,1)) || ~isempty(find(dz>dzMax,1)))
    % deal with out-of-bounds issue
    indMin=find(dz<dzMin);
    indMax=find(dz>dzMax);
    
    dz(indMin)=dzMin;
    dz(indMax)=dzMax;

    warning('palantir:SM_psf_lsqgauss:zOutOfRange','SM_psf_lsq1gauss encountered emission depths xyzEm(:,3) values outside of the valid model range.')
end

% PSF std rescaled to simulated wave-length
sig=lambda1/lambda0*polyval(Psigma,dz);

% generate positions at which photons are detected
xyDet = xyzEm(:,1:2)+(sig*[1 1]).*randn(Nph,2);
