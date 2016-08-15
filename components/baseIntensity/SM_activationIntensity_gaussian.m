% I=SM_activationIntensity_gaussian(opt)
%
% generate a basal intensity (photons/time unit) for a single fluorophore,
% based on the fields in the opt structure, which is the .baseIntensity
% field in the runinput file. Thismodel returns the same intensity
% opt.itensity every time.
%
% input     : option struct opt, with fields 
%             opt.intensity_mean
%             opt.intensity_std
% output    : I is a Gaussian with mean opt.mean and standard
%             deviation opt.std. Negative intensities are
%             discarded and resampled, so this model is only good if the
%             mean is considerably larger than the standard deviation.
%
% ML 2014-12-17

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_activationIntensity_gaussian.m, part of the SMeagol package
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
function I=SM_activationIntensity_gaussian(opt)
if(nargin==0) % return default options struct
    opt=struct('mean',1e5,'std',1e3);
    I=opt;
    return;
end

for n=1:10000
    I=opt.mean+opt.std*randn;
    if(I>=0)
        return
    end
end

error('SM_activationIntensity_gaussian: could not generate non-negative intensity in 10000 trials. Choose better parameters.')


