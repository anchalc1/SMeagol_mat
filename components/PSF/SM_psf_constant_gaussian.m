% function [xyDet,sigma]=SM_psf_lsq1gauss(t,xyzEm,opt)
%
% Simple PSF model for in-plane diffusion, it just adds a zero-mean
% Gaussian noise with a specified standard deviation to every
% x,y-coordinate.
%
% t         - emission times (one per photon, not used here)
% xyzEM     - positions of photon emitters (one per photon), one
%             xyz-triplet per row.
% xyDet     - positions of detected photons on the camera chip, randomly
%             generated using the PSF model.
%
% opt       - PSF options
% opt.sigma - standard deviation of the Gaussian PSF
%
% The options object opt here corresponds to the opt.psf structure defined
% in the Palantir runinput file. Calling the function without input
% arguments returns a default options struct with sigma=130;
%
% M.L. 2015-01-15

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_constant_gaussian, simulated PSF in the SMeagol package
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
function [xyDet,sig]=SM_psf_constant_gaussian(~,xyzEm,opt)

if(nargin==0) % return default options struct
    opt=struct('sigma',130);
    xyDet=opt;
    sig=[];
    return;
end
sig=opt.sigma;
Nph=size(xyzEm,1);
% generate positions at which photons are detected
xyDet = xyzEm(:,1:2)+sig*randn(Nph,2);
