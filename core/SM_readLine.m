% [tline,id,s,x]=SM_readLine3(fid,A,b,tScale,speciesStr,lastSpeciesIsDegraded,ignoreThis)
%
% read one line from an input file (trajectory or reactions), and extract
% time t, particle id numbers id, states s, and positions x (as an N*3
% array). Empty vectors are returned if no line is left to read.
%
% fid : file opinter to read from (from fopen('...','r')
% A,b ; transformation matrix and offset vector for conversion from voxel
%       to actual coordinates. See opt.trj.A, opt.trj.b
% tScale: conversion factor from input time (read from data) to output time
%         (tline).
% speciesStr : cell vector of species names, in the same order as diffusion
%              constants
% delastSpeciesIsDegraded : if true, the last entry in speciesString is
%                           interpreted as the degraded state s=-1.
%                           Default=false (no special treatment of a
%                           degraded state).
% ignoreThis : boolean indicator vector. A particle is omitted from output
%              if its pid-number id satisfies ignoreThis(id)=true.
%              Particles not mentioned in the trackThis vector, e.g., if
%              id > length(trackThis), are not omitted. 
%              (Default: no omissions)
%                           
%
% ML 2013-11-26

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_readLine.m, read single line of input data for the SMeagol package.
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
function [tline,id,s,x]=SM_readLine(fid,A,b,tScale,speciesStr,lastSpeciesIsDegraded,ignoreThis)

if(~exist('lastSpeciesIsDegraded','var') || isempty(lastSpeciesIsDegraded))
    lastSpeciesIsDegraded=false;
end
degradedIndex=0;
if(lastSpeciesIsDegraded)
    degradedIndex=length(speciesStr);
end
if(~exist('ignoreThis','var'))
    ignoreThis=[];
end

%%%try % read the next trajectory entry
linestr=fgetl(fid);
if(linestr~=-1) % foo == -1 means eof reached
    spcind=[strfind(linestr,' ') length(linestr)+1]; % indices to spaces btw entries, including eol
    
    % remove double spaces
    ind2spc=spcind(diff(spcind)==1);
    if(~isempty(ind2spc))
       linestr=linestr(setdiff(1:length(linestr),ind2spc)); 
       spcind=[strfind(linestr,' ') length(linestr)+1];
    end
    
    Nparticles= (length(spcind)-1)/5; % number of particle entries on this line
    tline=tScale*str2double(linestr(1:(spcind(1)-1))); % time-stamp
    
    id=[];s=[];x=[];
    % extract particle information    
    for nid=1:5:(5*Nparticles) % spaces just before all particle ids in linestr
        iid0=spcind(nid)+1;
        iid1=spcind(nid+1)-1;
        idNew=str2double(linestr(iid0:iid1));
        
        % identify species
        sid0=spcind(nid+1)+1;
        sid1=spcind(nid+2)-1;
        sNew=find(strcmp(linestr(sid0:sid1),speciesStr),1);
        if(lastSpeciesIsDegraded && sNew==degradedIndex)
            sNew=-1;
        end
        if(isempty(sNew))
            warning(['unidentified species ' linestr(sid0:sid1) ', ignoring.'])
            continue
        end
        if( idNew<=length(ignoreThis) && ignoreThis(idNew) ) % ignore particles if told so
            continue
        end
        id(end+1)=idNew;
        s(end+1)=sNew;
            
        % extract coordinates
        xid0=spcind(nid+(2:4))+1;
        xid1=spcind(nid+(3:5))-1;
        x(1:3,end+1)=[str2double(linestr(xid0(1):xid1(1))); ...
            str2double(linestr(xid0(2):xid1(2))); ...
            str2double(linestr(xid0(3):xid1(3)))];
    end
    if(~isempty(x))
        x=(A*x+b*ones(1,size(x,2)))'; % convert to actual positions
    end
else
    tline=[];
    id=[];
    s=[];
    x=[];
    return
end
