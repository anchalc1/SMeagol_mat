% xyDet=SM_psf_noPSF(t,xyzEm,opt)
%
% Dummy PSF function that can be used when movies are not of interest; it
% simply returns the input coordinates.
%
% input:
% t         - photon emision times (one per photon, not used here)
% xyzEM     - positions of photon emitters, one xyz triplet per row for
%             each emitted per photon.
% opt       - PSF options (not used here, but something must be passed).

% output:
% xyDet     - positions of detected photons on the camera chip, randomly
%             generated using the PSF model.
%
%
% M.L. 2015-01-15

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_psf_noPSF, simulated (dummy) PSF in the Smeagol package
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
function xyDet=SM_psf_noPSF(~,xyzEm,~)

if(nargin==0) % return default options struct
    opt=struct;
    xyDet=opt;
    return;
else
    xyDet=xyzEm(:,1:2);
end
