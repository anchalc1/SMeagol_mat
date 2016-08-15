% photophys=SM_fluo_multistate_blinking(opt)
% A phophysical model with multiple dark states.
%
% input: options struct, with fields
% opt.tau       : mean dwell times of the dark states
% opt.weights   : relative weights (>=0) of the times tau in the multistate
%                 dark time distribution.
% opt.bleachTime  : inverse bleaching rate from all states when
%                   illuminated.
% opt.fluoFraction: fraction of time spent in the fluorescent state (this
%                   sets the overall conversion rate to the dark states).
%                   (More precisely, this number refers to the steady state
%                   occupancy in the illuminated states with no bleaching.)
% 
% Output: Palantir photophysics opt struct, with fields 
% opt.emissionFactor
% opt.kb      = {kb1, kb2}                  Bleaching rates
% opt.Q       = {Q1 Q2};                    Internal state dynamics
%                                           {with,without} excitation laser.
%
% With no input, the function returns some default options, which have a
% power-law like 1/t^2 dwell time density in the range [1e-5 1e-2] s.

function photophys=SM_fluo_multistate_blinking(opt)

if(nargin==0) % return default options struct
    % a simple photophysics model with a single bright state and no
    % bleaching
    opt=struct;
    opt.tau    =[3.73e-06 1e-05 2.68e-05 7.20e-05 1.93e-4 5.18e-4 1.39e-3 0003.73e-3];  % blinking dwell times
    opt.weights=[0.636    0.228 0.0843   0.0328   0.0111  0.005   0.00115 0.00106   ];    % relative blink time weights
    opt.fluoFraction=0.75;       % fraction of time spent in fluorescent state
    opt.bleachTime=0.1;         % overall bleaching time when illuminated.
    photophys=opt;
    return;
end

% if an option was given, compute the kinetics parameters
tauMean=sum(opt.weights.*opt.tau)/sum(opt.weights);          % mean blinking time
bTotal =(1-opt.fluoFraction)/opt.fluoFraction/tauMean;  % overall blink rate

N=1+length(opt.tau); % number of states

photophys=struct;
photophys.emissionFactor=[1 zeros(1,N-1)];
photophys.kb={ones(1,N)/opt.bleachTime,zeros(1,N)}; % no bleaching
photophys.Q=cell(1,2);

Q=zeros(N,N);
Q(2:N,1)=reshape(1./opt.tau,N-1,1);   % recovery from blinking
photophys.Q{2}=Q;                     % without illumination: only recovery
Q(1,2:N)=bTotal*opt.weights;          % blinking rates
photophys.Q{1}=Q;                     % with illumination: recovery&blinking



