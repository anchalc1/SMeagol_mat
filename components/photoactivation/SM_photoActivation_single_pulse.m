% tActivation=SM_photoActivation_single_pulse(t0,activation,sample)
% 
% Simulate a photoactivation event history starting in an inactivated state
% at time t0, given a SM options object SMopt. The output is an activation
% time >= t0.
%
% SM_photoActivation_single_pulse models photo-activation by one
% photoactivation event at time tP, where existing fluorophores are
% activated to state sStart=1 with probability Pa, and particles created
% after that are not activated at all. 
%
%
% input:
% t0    : time at which to start the activation simulation.
% activation : SMeagol activation options (opt.activation). Required fields:
%   ta  : time for activation event
%   Pa  : activation probability at the activation event. So, if t0<=ta,
%         then tActivation=ta with probability Pa, and Inf otherwise. If
%         t0>ta, then tActivation=inf (no activation).
% sample  : Standard SMeagol sampling options. Not usewd here.
%
% output:
% tActivation   : photoactivation time
%
% ML 2015-12-11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_photoActivation_single_pulse.m, simulate photoactivation time as part
% of  the SMeagol package
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
function tActivation=SM_photoActivation_single_pulse(t0,activation,sample)
% return default options if no input arguments were given

if(nargin==0) % return default options struct
    activation=struct('ta',0,'Pa',0.01);
    tActivation=activation;
    return;
end

% activation 
if(t0<=activation.ta && rand<=activation.Pa)
    tActivation=activation.ta;
else
    tActivation=Inf;
end
