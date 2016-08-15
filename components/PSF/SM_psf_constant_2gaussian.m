% function [xyDet,sigma]=SM_psf_constant_2gaussian(t,xyzEm,opt)
%
% Slightly less simple PSF model for in-plane diffusion that models the PSF
% as a superposition of two Gaussians with the same mean and no
% z-dependence.
%
% t         - photon emision times (one per photon, not used here)
% xyzEM     - positions of photon emitters (one per photon).
% xyDet     - positions of detected photons on the camera chip, randomly
%             generated using the PSF model.
%
% opt       - PSF options
% opt.sigma - standard deviations of the Gaussian PSFs
% opt.frac1 - relative weight of the Gaussian corresponding to sigma(1),
%             should be between 0 and 1.
%
% The options object opt here corresponds to the opt.psf structure defined
% in the SMeagol runinput file. Calling the function without input
% arguments returns a default options struct with 
% sigma=[130 395] nm, and frac1=0.65
%
% The two-Gaussian approach was inspired by 
% Zhu, L., Zhang, W., Elnatan, D. & Huang, B. 
% Faster STORM using compressed sensing. Nat Meth 9, 721–723 (2012).
%
% M.L. 2015-01-15

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_constant_2gaussian, simulated PSF in the SMeagol package
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
function [xyDet,sig]=SM_psf_constant_2gaussian(~,xyzEm,opt)

if(nargin==0) % return default options struct
    opt=struct('sigma',[130 395],'frac1',0.65);
    xyDet=opt;
    sig=[];
    return;
end
Nph=size(xyzEm,1);
sig=opt.sigma(1)*ones(Nph,1);
sig(rand(Nph,1)>opt.frac1)=opt.sigma(2);

% generate positions at which photons are detected
xyDet = xyzEm(:,1:2)+(sig*ones(1,2)).*randn(Nph,2);
