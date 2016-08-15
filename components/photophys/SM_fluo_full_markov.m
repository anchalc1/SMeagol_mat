% photophys=SM_fluo_full_markov(opt)
%
% A generator of SMeagol photophysics models from input parameters opt,
% which could be the opt.photophys field of the Palantir runinput options.
%
% In this case, it is a dummy generator that simply returns the appropriate
% fields of the input struct, which is cumbersome to handle explicitly, but
% useful as a template for creatig your own blink&bleach models.
% 
% The input (and output) fields are 
% opt.emissionFactor    relative emission intensity of the illuminated
%                       states. 1*N row-vectors, where N is the number of
%                       photophysical states.
% opt.kb   = {kb1 kb2}; bleaching rates of the states when illuminated (kb1) 
%                       and not illuminated (kb2). N*1 column-vectors,
%                       where N is the number of photophysical states.
% opt.Q    = {Q1   Q2}; Interconversion rate matrices (continuous time)
%                       when illuminated (Q1) and not illuminated (Q2). N*N
%                       matrices. Only off-diagonal elements are used in
%                       the simulations, and normalizations are enforced
%                       elsewhere: Q(i,i)=-sum( Q(i,[1:i-1 i+1:end]);
%
% With no input, the function returns default parameters describing an
% ideal 1-state fluorophore that does not blink or bleach.
%
% To program your own models, make them return the same output fields as
% this function.
%
% M.L. 2014-12-18

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_full_photophysics, template for photophysics models in SMeagol
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
function photophys=SM_full_photophysics(opt)

if(nargin==0) % return default options struct
    % a simple photophysics model with a single bright state and no
    % bleaching
    opt=struct;
    opt.emissionFactor=1;
    opt.kb={[0]  [0]};
    opt.Q={[0]  [0]};
    photophys=opt;
    return;
end
% if an option was given, just pass the photophysics fields along:
photophys=struct;  % create empty output structure
photophys.emissionFactor=opt.emissionFactor;
photophys.kb=opt.kb;
photophys.Q=opt.Q;
