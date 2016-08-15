% A gamma-distribution model for the EM multiplication process, using the
% large-gain/high signal approximation (eqs. 13-14) of Hirsch et al [1].
%
% usage: 
% nOut = SM_EMgain_hirsch(nIn,gain)
% 
% nOut  : number of output electrons from the simulated EM gain process
% nIn   : number of input electrons
% gain  : EM gain, here defined as the average number of output electrons
%         per input electron.
%
% Model: nOut = nIn - 1 + xi, where xi ~ Gamma-distributed with shape
% parameters nIn and shape parameter theta= gain-1+1/nIn, rounded to
% nearest integer. With this definition, we have
% 
% <nOut> = gain * nIn, Var(nOut) = nIn * gain^2, i.e., an excess noise
% factor of exactly sqrt(2).
%
% 1. Hirsch, M., R.J. Wareham, M.L. Martin-Fernandez, M.P. Hobson, and D.J.
% Rolfe. 2013. A Stochastic Model for Electron Multiplication
% Charge-Coupled Devices – From Theory to Practice. PLoS ONE. 8: e53671.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_EMgain_hirsch.m, part of the SMeagol package
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
function nOut = SM_EMgain_hirsch(nIn,gain)
theta=gain+1./nIn-1;
xi=round(gamrnd(nIn,theta));
nOut=nIn-1+xi;
