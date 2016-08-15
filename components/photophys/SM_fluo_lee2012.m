% photophys=SM_full_photophysics(opt)
%
% A generator of Palantir photophysics models from input parameters opt,
% which could be the opt.photophys field of the Palantir runinput options.
%
% In this case, the fluorophore model is one with Markovian blinking and
% bleaching, inspired by the mEos2 and Dendra2 models presented by Lee et
% al (2012,*), but with the added assumption that blinking and bleaching
% only takes place when illuminated by an excitation laser (also note that
% some of these rates presumably depend on the applied laser instensities). 
%
% The activation kinetics is not part of this model, but in the separate
% activation model.
%
% States: 1) fluorescent, 2,3) dark, and -1) bleached. 
% Illuminated : blink, recover, bleach
%     ------ kr1 ----->
% (2) <----- kd*a ----- (1) 
%                       (1) --- kb ---> (-1)
% (3) ------ kr2 -----> (1)
%     <--- (1-a)*kd ---
%
% Non-illuminated: only recovery    
% (2) ------ kr1 -----> (1) 
%                       (1)             (-1)
% (3) ------ kr2 -----> (1)
%
% The rates are parameterized in terms of the mean lifetimes t1,t2,t3 of
% states 1-3 (under illumination), the bleaching probability pB=kb/(kb+kd),
% and the conditional dark state splitting probability a. In this notation,
% the kinetics reported in [*] would be  as follows:
% mEos2   : t123=[0.076 2.5  0.064], pB=0.41, a=0.81
% Dendra2 : t123=[0.051 0.63 0.056], pB=0.82, a=0.24
%
% 
% Input: opt struct, with fields
% opt.t123  : mean life-times of states 1-3 under illumination. 
%             (All mean life-times > 0).
% opt.pB    : conditional bleaching probability after an illuminated
%             fluorescent state. (0 <= opt.pB <= 1).
% opt.a     : conditional dark state probability (to end up in state 2),
%             (0 <= opt.a <= 1). 
%
% Output: Palantir photophysics opt struct, with fields 
% opt.emissionFactor = [1 0 0]; (States 2,3 are dark).
% opt.kb      = {[pB/t123(1);0;0] [0;0;0]}; bleaching only from illuminated
%                                           fluorescent state. 
% opt.Q       = {Q1 Q2};                    Internal state dynamics
%                                           {with,without} excitation laser.
%
% With no input, the function returns default options, corresponding to the
% mEos2-inspired parameters, in units of seconds.
%
% M.L. 2014-12-19
%
% (*) Lee, S.-H., Shin, J. Y., Lee, A. & Bustamante, C. Counting single
% photoactivatable fluorescent molecules by photoactivated localization
% microscopy (PALM). PNAS 109, 17436–17441 (2012). 
% www.pnas.org/cgi/doi/10.1073/pnas.1215175109

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_fluo_lee2012, a photophysics model in the SMeagol package
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

%% actual code
function photophys=SM_fluo_lee2012(opt)

if(nargin==0) % return default options struct
    % a simple photophysics model with a single bright state and no
    % bleaching
    opt=struct('t123',[0.076 2.5  0.064],'pB',0.41,'a',0.81);
    photophys=opt;
    return;
end
% if an option was given, compute the kinetics parameters
kB=opt.pB/opt.t123(1);
kD=(1-opt.pB)/opt.t123(1);
kr1=1/opt.t123(2);
kr2=1/opt.t123(3);
a=opt.a;

photophys=struct;  % create empty output structure
photophys.emissionFactor=[1 0 0];
photophys.kb={[kB 0 0]',[0 0 0]'};
Q1=[0    kr1 kr2;
    a*kD   0   0;
(1-a)*kD   0   0];
Q2=[0   kr1 kr2;
    0     0   0;
    0     0   0];
photophys.Q={Q1,Q2};
