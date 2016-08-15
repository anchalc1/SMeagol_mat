% function I=SM_activationIntensity_uniform(opt)
%
% generate a basal intensity (photons/time unit) for a single fluorophore,
% based on the fields in the opt structure, which is the .baseIntensity field
% in the runinput file. Thismodel returns the same intensity opt.itensity
% every time.
%
% input     : option struct opt, with field opt.intensity
% output    : I=opt.intensity. If no input parameters are given, a default
%             opt struct is returned, with intensity=1e5;
%
% ML 2014-12-17

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_activationIntensity_uniform, in the SMeagol package.
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

function I=SM_activationIntensity_uniform(opt)
if(nargin==0) % return default options struct
    opt=struct('intensity',1e5);
    I=opt;
    return;
end

I=opt.intensity;
