% tActivation=SM_photoActivation_instant_numbered(t0,activation,sample)
% 
% Simulate a single photoactivation event history starting in an inactivated
% state at time t0, given some activation and sample options. The output is
% an activation time >= t0. This particular model activates a predetermined
% set of particles (defined through the order in which they are called) to
% tActivation=t0, while other particles are not activated. This is good
% for, e.g., systematically looping through all fluorophores in a big
% simulation.
%
% CAUTION: particle counter is a persistent variable (see matlab help text
% on peristent), and thus it does not reset automatically when a SMeagol
% simulation starts. Calling this function without input arguments will
% reset it, but this must be done manually, or by opening the
% photoactivation parameters dialog in the SMeagol gui.
%
% input:
% t0          : time at which the activation process is started. 
% activation  : SMeagol activation options (opt.activation). 
% activation.pnum = [n1 n2 .. nN]; which particles to activate. Default
%               [1], i.e., activate the first particle.
% sample      : SMeagol sampling options (opt.sample). No fields used.
%
% output:
% tActivation : photoactvation time
%
% ML 2016-03-08

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_photoActivation_instant_numbered, part of the SMeagol package
% ========================================================================= 
% Copyright (C) 2016 Martin Lindén and Johan Elf
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
function tActivation=SM_photoActivation_instant_numbered(t0,activation,~)

persistent pcounter

% return default options if no input arguments were given
if(nargin==0) % return default options struct
    activation=struct();
    activation.pnum=1;
    tActivation=activation;
    pcounter=0;
    disp('SM_photoActivation_instant_numbered resetting particle counter!')
    return;
end

if(isempty(pcounter))
    pcounter=0;
    disp('SM_photoActivation_instant_numbered resetting particle counter!')
end
pcounter=pcounter+1;
%disp([' SM_..._numbered: ' int2str(pcounter) ' at ' num2str(t0)])

if(ismember(pcounter,activation.pnum)) 
    tActivation=t0;
else
    tActivation=Inf;
end
