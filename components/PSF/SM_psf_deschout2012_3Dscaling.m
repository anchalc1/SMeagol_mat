% [xyDet,sigma]=SM_psf_deschout2012_3Dscaling(t,xyzEm,opt)
%
% Explicit PSF using the Gaussian model of Deschouet al, 
% Deschout, H., Neyts, K. & Braeckmans, K. The influence of movement on the
% localization precision of sub-resolution particles in fluorescence
% microscopy. Journal of Biophotonics 5, 97–109 (2012).  
% 
% The model is Gaussian, with a z-dependent width/std s(z), given by
% s(z)^2 = s0^2 + (z-zp)^2*(s1/s0)^2, where
% s0 = with at focal plane (probably proportional to the wavelength)
% zp = position of focal plane (defined as zero here)
% s1 = a characteristic scaling length, given by
% s1 = opt.lambda/(4*pi*opt.n)
%
% output: 
% xyDet     - positions of detected photons on the camera chip, randomly
%             generated using the above PSF model.
% sigma     - list of PSF std used to generate the corresponding positions
%
% input:
% t         - emission time for each photon (not used here).
% xyzEM     - positions of photon emitters (one per photon).
% opt       - PSF options
% opt.s0    - PSF standard deviation in the focal plane.
% opt.lambda- emitted light wavelength
% opt.n     - index of refraction in the sample
%
% The options object opt here corresponds to the opt.psf structure defined
% in the Palantir runinput file. Calling the function with no input
% parameters returns a default options struct.
%
% M.L. 2014-12-16

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_deschout2012_3D, simulated PSF in the SMeagol package
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
function [xyDet,sigma]=SM_psf_deschout2012_3Dscaling(~,xyzEm,opt)

 % return default options struct if no input given
if(nargin==0)
    opt=struct('s0',130,'lambda',560,'n',1.33);
    xyDet=opt;
    sig=[];
    return;
end

s0= opt.s0;   % static psf std in focal plane
z0= 4*pi*opt.n/opt.lambda*s0^2; % characterisitc length for out-of-focus psf widening
                        % s(z) = s0*sqrt(1+ (z-zp)^2/z0^2), 
                        % from Deschout et al 2012.
zFP=0;    % reference height of the focal plane

% emitted photons
Nph=size(xyzEm,1); % number of emitted photons

% positions at which photons are detected
xyDet = xyzEm(:,1:2)+(s0*sqrt(1+((xyzEm(:,3)-zFP)/z0).^2)*[1 1]).*randn(Nph,2);

sigma=s0*sqrt(1+((xyzEm(:,3)-zFP)/z0).^2); % compute the std too, if neededs
