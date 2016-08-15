% tActivation=SM_photoActivation_Pa_basal(t0,activation,sample)
% 
% Simulate a single photoactivation event history starting in an inactivated
% state at time t0, given a SM options object SMopt. The output is an
% activation time > t0.
%
% SM_photoActivation_Pa_basal models photo-activation by 1) the activation
% laser, and 2) spontaneous activation, with different rates when the
% illumination laser is on or off. The activated fluorophore starts in
% state sStart=1.
%
% Details:
% 1) A constant activation probability activation.Pa every time the
% activation laser is fired, which occurs every activation.ta. An
% additional random time with exponential distribution and mean value 
% activation.td (can be zero) is added, to enable additional randomness to
% the simulation.
%
% 2) Spontaneous activation with an activation rate activation.ka(I) per
% unit time, with I=1 when illuminated, and I=2 otherwise.
%
% input:
% t0    : time at which the activation process is started. 
% activation : SMeagol activation options (opt.activation). Required fields:
%   t1  : time for first activation laser pulse
%   ta  : time inteval for activation laser, with first pulse at t=t1.
%         ta=inf is an error, the special case of a single activation pulse
%         can be set by setting t1 = pulse time, ta > total simulation
%         time (but finite).
%   Pa  : activation probability every time the activation laser is fired.
%   td  : average delay time for activation  with activation laser.
%         td=0 gives no delay.
% ka=[ka1 ka2] : spontaneous activation rates [with without] illumination
%         laser. 
% sample  : Standard SMeagol sampling options. Used fields:
%           dt  : length of sampling time interval, starting at t=0.
%           tE  : exposure/illumination/aquisition interval. 0 < tE <= dt
%           is required. During each aquisition interval, ka(1) is the
%           activation rate for 0 < t < tE, and the ka(2) is used during tE
%           < t < dt.
%
% output:
% tActivation   : photoactivation time
%
% ML 2014-12-17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_photoActivation_Pa_basal.m, simulate photoactivation time as part of 
% the SMeagol package
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
function tActivation=SM_photoActivation_Pa_basal(t0,activation,sample)
% return default options if no input arguments were given

if(nargin==0) % return default options struct
    activation=struct('t1',0,'ta',50e-3,'Pa',0.01,'td',0,'ka',[0 0]);
    tActivation=activation;
    return;
end

%% activation by activation laser, probability activation.Pa per flash
na0=ceil((t0-activation.t1)/activation.ta); % number of activation pulses before t0
na0=max(0,na0); % no pulses before activation.t1
ta1=Inf;
if(activation.Pa>0)
    if(activation.Pa<1)
        na1=ceil(log(1-rand)/log(1-activation.Pa)); % number of activation attempts needed for success
    elseif(activation.Pa>=1)
        na1=1; % instant success 
    end
    
    ta1=activation.t1+(na0+na1-1)*activation.ta; % time of the activating flash
    ta1=ta1-activation.td*log(1-rand);           % add delay time
end
clear na0 na1
        
if(sum(activation.ka)==0) % no spontaneous activation. We are done.
   tActivation=ta1;
   return
end
%% consider activation by spontaneuos activation
t1=sample.tE;         % illumination time per cycle
t2=sample.dt-sample.tE; % unilluminated time
dt=t1+t2;
ka=activation.ka;
kMean=(t1*ka(1)+t2*ka(2))/dt; %mean activation rate

u=rand;
na_min=floor( -log(1-u)/dt/kMean); % number of unsuccessful activation periods

% need to figure out when activation happens during the successful period
tt=[0 t1 dt dt+t1 2*dt];
kaFrac=cumsum([0 t1*ka(1) t2*ka(2) t1*ka(1) t2*ka(2)])/dt/kMean;

ur=-log(1-u)/dt/kMean-na_min; % integral reminder
ta0=mod(t0,dt); % where in the sampling period does t0 start?
ua0=interp1(tt,kaFrac,ta0);

%taa=interp1(kaFrac,tt,ua0+ur); % activation time, after the start of the last activation period
% this does not work if ka has some zero elements. Instead, we interpolate
% manually, and consistently choose the highest activation time.
i0=find(kaFrac<=ua0+ur,1,'last');
taa=tt(i0)+(tt(i0+1)-tt(i0))*(ua0+ur-kaFrac(i0))/(kaFrac(i0+1)-kaFrac(i0));

ta2=t0+dt*na_min-ta0+taa;

% see which type of activation happened first.
tActivation=min([ta1 ta2]);
