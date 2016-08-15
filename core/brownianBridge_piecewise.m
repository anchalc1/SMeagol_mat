% xi=brownianBridge_piecewise(t,x,D,ti)
% 
% Performs Brownian bridge target to times ti, given positions x at 
% times t, and diffusion constants D. Index conventions:
% x(m,:) is position at time t(m), and D(m) is the diffusion constant from
% t(m) < t <= t(m+1). ( D(end) is not used). 
%
% Both t, and ti, must be sorted in ascending order, and ti must be inside
% the interval spanned by t, 
% t(1) <= ti(n) <= t(end).
%
% ML 2013-11-15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% brownianBridge_piecewise.m, part of the SMeagol package.
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


function xi=brownianBridge_piecewise(t,x,D,ti)

if(length(t)<2)
    error('brownianBridge_piecewise: need an interpolation interval')   
end
if(~issorted(ti))
    error('brownianBridge_piecewise: target times ti must be in ascending order')   
end
if(~issorted(t))
    error('brownianBridge_piecewise: boundary times t must be in ascending order')   
end
if(ti(1) < t(1) || ti(end)>t(end))
    error('brownianBridge_piecewise: cannot interpolate outside time interval')
end

if(sum(diff(t)<0)>0)
   error('interpolation times must be increasing') 
end

dim=size(x,2);


% start of first Brownian bridge
iStart=find(t>ti(1),1);

xi=zeros(length(ti),dim);
ii=1; % interpolation point index
iiMax=length(ti);
for iNext=iStart:length(t)
    % next known position
    xNext=x(iNext,:);
    tNext=t(iNext);
    
    % last known position, time, and diffusion constant
    xLast=x(iNext-1,:);
    tLast=t(iNext-1);
    Dcurrent=D(iNext-1);
    
    while(ii<=iiMax && ti(ii)<=tNext)        
        dtNormalized=(ti(ii)-tLast)/(tNext-tLast);
        xi(ii,:)=xLast+(xNext-xLast)*dtNormalized... % linear interpolation
            +sqrt(2*Dcurrent*(ti(ii)-tLast)*(1-dtNormalized))*randn(1,dim); % random part
        
        % update last known position and time
        xLast=xi(ii,:);
        tLast=ti(ii);
        
        ii=ii+1;
    end    
end
