% photophys=SM_fluo_only_bleach(opt)
%
% A generator of Palantir photophysics models from input parameters opt,
% which could be the opt.photophys field of the Palantir runinput options.
%
% In this case, the fluorophore model is one with simple bleaching dynamics
% (no blinking), and the user simpliy specifies the average lifetime of the
% fluorescent state. Bleaching supposedly only takes place during
% illumination, and so this can be used to simulate a simple trade-off
% between long exposure and many images.
% 
% Input: opt struct, with fields
% opt.bleach_time   : mean life-time tB of fluorescent state when the
%                     fluorophore is illuminated. Must be
%                     positive. Default=inf (no bleaching).
% Output: Palantir photophysics opt struct, with fields 
% opt.emissionFactor = 1; (Intensity as specified by the activation model).
% opt.kb      = {1/tB 0}; bleaching only when illuminated.
% opt.Q          = {0 0}; No internal dynamics except bleaching.
%
% With no input, the function returns default parameters describing an
% ideal 1-state fluorophore that does not blink or bleach.
%
% M.L. 2014-12-18

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_full_photophysics, a photophysics model in SMeagol
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
function photophys=SM_fluo_only_bleach(opt)

if(nargin==0) % return default options struct
    % a simple photophysics model with a single bright state and no
    % bleaching
    opt=struct('bleach_time',inf);
    photophys=opt;
    return;
end
% if an option was given, just pass the photophysics fields along:
tB=opt.bleach_time;
if(tB>0)
    photophys=struct;  % create empty output structure
    photophys.emissionFactor=[1];
    photophys.kb={1/tB,0};
    photophys.Q={0,0};
else
   error('SM_fluo_only_bleach: opt.bleach_time > 0 is required.')
end
