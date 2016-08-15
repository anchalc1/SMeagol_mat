% testresult=SM_mesoRD_tracking_test(mesoRDstring)
% a simple test of mesoRD tracking capability: run a simple SBML model and
% see if tracking files are produced.
%
% input:
% mesoRDstring : system-specific command to start mesord. 
%                Default: look in the expected place.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_mesoRD_tracking_test.m, part of the SMeagol package.
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
function testresult=SM_mesoRD_tracking_test(mesoRDstring)

% folders
startFolder=pwd;
[coreFolder,~,~]=fileparts(mfilename('fullpath'));
% the place where MesoRD should be:
if(ispc)
  mesoRDstring=fullfile(coreFolder,'.','..','mesord','mesord.exe');
else
  mesoRDstring=fullfile(coreFolder,'.','..','mesord','mesord');
end
if(exist(mesoRDstring,'file')~=2)
    warning(['Cannot find MesoRD in the correct place : ' mesoRDstring])
    testresult=false;
    return
end
%% test the tracking capabilities
outputfiles={'geometry.pov','geometry.txt','species_A.txt',...
    'species_B.txt','molecule_tracking_log.txt','time.txt',...
    'reactions.txt','trajectories.txt'};

% remove existing output files
for k=1:length(outputfiles)
    if(exist(fullfile(coreFolder,outputfiles{k}),'file'))
        delete(fullfile(coreFolder,outputfiles{k}));
    end
end
cd(coreFolder)
% hack to work around the use of pathnames with spaces: add " around the
% path to mesord
system(['"' mesoRDstring '" -i 1 -I 50 -c 1 -C -1 -E -g SM_mesoRD_tracking_test.xml -t 0.1 -q 0.20 um -K -x 0.1 -w 0.05']);

cd(startFolder)
fprintf('\n')
%rm geometry.* reactions.txt trajectories.txt time.txt species_* molecule_tracking_log.txt

%%% check presence of tracking output:
disp('-----------------------------')
if(exist(fullfile(coreFolder,'molecule_tracking_log.txt'),'file'))
    testresult=true;
    disp([' MesoRD tracking works using : ' mesoRDstring])
    
else
    testresult=false;
    disp([' MesoRD tracking does not work using : ' mesoRDstring])
    disp(' Try another location, and make sure that the MOLECULE_TRACKING option is enabled.')
end
disp('-----------------------------')
% remove existing output files again
for k=1:length(outputfiles)
    if(exist(fullfile(coreFolder,outputfiles{k}),'file'))
        delete(fullfile(coreFolder,outputfiles{k}));
    end
end
