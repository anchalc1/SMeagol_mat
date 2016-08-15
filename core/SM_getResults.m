% [res,opt]=SM_getResults(runinput)
%
% Retrieved options struct and loads the results mat file if it is available
% (if not, a warning is issued, and an empty matrix).
% The optional runinput can be a file name (with path) or an options
% struct. If not given, a file dialog for a runinput file is opened.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_getResults.m, read simulation output in the mesoSM package.
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining
%  it with Matlab or any Matlab toolbox, the licensors of this Program
%  grant you additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code

function [res,opt]=SM_getResults(runinput)

if(exist('runinput','var'))
    opt=SM_getOptions(runinput);
else
    opt=SM_getOptions();
end

resultFile=fullfile(opt.runinputroot,opt.output.resultFile);

if(exist(resultFile,'file'))
    res=load(resultFile);
else
    res=[];
    warning(['No results found for ' opt.runinputfile ' (yet).'])
end
