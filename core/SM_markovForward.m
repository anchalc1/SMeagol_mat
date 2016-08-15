% [s1,tsHistory]=SM_markovForward(Q,s0,t0,t1,QT)
%
% Advance a markov chain with rate matrix Q from t=t0 to t=t1, starting in
% state s0, and report all transitions that occur. Initial and final states
% are also included in tsHistory.
% The rate convention is q(i,j) = r{ i -> j }, that is, dp/dt = p*Q, where
% p is a row matrix. Row normalization, sum(Q(k,:)) = 0, is enforced by
% ignoring the diagonal elements of Q.
%
% QT: optional termination rates (e.g., bleaching), which governs
% irreversible transitions to state s=-1. Should be N*1 vector, or scalar
% (same rate for all states). If omitted, no terminations occur.
%
% Martin Lindén 2013-11-27 (updated to s=-1 for bleaching)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_markovForward, simulate Markov process in continuous time for the 
% SMeagol package.
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


function [s1,tsHistory]=SM_markovForward2(Q,s0,t0,t1,QT)

N=size(Q,1); % number of states
sT=N+1; % termination state
if(exist('QT','var'))
    [rT,cT]=size(QT);
    if(length(QT)<N)
        QT=QT*ones(N,1);
    else
        
    end
else
    QT=zeros(N,1);
end

Qtot=[Q.*(1-eye(size(Q))) QT;zeros(1,N+1)];
    
sCurrent=s0;
tCurrent=t0;
tsHistory=[t0 s0];
while tCurrent < t1
   % assemple all transition rates   
   rOut=sum(Qtot(sCurrent,:));     % total out-rate
   pOut=cumsum(Qtot(sCurrent,:)/rOut);     % jump probabilities
   if(rOut>0)
       tNext=tCurrent-1/sum(rOut)*log(1-rand); % time until next event
       sNext=find(rand<pOut,1);
       if(tNext <t1)
           tCurrent=tNext;
           sCurrent=sNext;
           tsHistory(end+1,:)=[tNext sNext];
       else
           tCurrent=t1;
       end
   else
      tCurrent=t1;      
   end
    
    
end
s1=sCurrent;
tsHistory(end+1,:)=[t1 s1];

% change termination state to zero
tsHistory(tsHistory(:,2)==sT,2)=-1;
s1=tsHistory(end,2);
