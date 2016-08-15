% tActivation=SM_photoActivation_Pa_instant(t0,activation,sample)
% 
% Simulate a single photoactivation event history starting in an inactivated
% state at time t0, given some activation and sample options. The output is
% an activation time >= t0. In this case, we use the simple model
% tActivation=t0, no free parameters.
%
% input:
% t0          : time at which the activation process is started. 
% activation  : SMeagol activation options (opt.activation). No required
%               fields. 
% sample      : SMeagol sampling options (opt.sample). No fields used.
%
% output:
% tActivation : photoactvation time
%
% ML 2014-12-17

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_photoActivation_instant, part of the SMeagol package
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
function tActivation=SM_photoActivation_instant(t0,~,~)
% return default options if no input arguments were given
if(nargin==0) % return default options struct
    activation=struct();
    tActivation=activation;
    return;
end
% the simplest photoactivation model of all!
tActivation=t0;
