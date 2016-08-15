function sysflag=SM_runMesoRD(mesoRDoptions,targetFolder)
% SM_runMesoRD(mesRDoptions,targetFolder)
% Call the mesord binary with a target folder and commandline options
% (including path to a model file) to see if tracking is supported.
% 
% input :
% mesoRDoptions : command-line options for mesoRD. '-h' will print the
%                 mesoRD help text, including a list of options. See also
%                 the mesoRD documentation.
% targetFolder  : folder from which mesoRD is called, which is also where
%                 any output is written. This is done using cd, so
%                 targetFolder can be either absolute or relative to the
%                 current path. Before ending, the current path is reset to
%                 the one from which SM_runMesoRD was called.
%                 default: pwd
% output : The output from mesoRD is written to the matlab command prompt.
% sysflag       : the output from the call to system() that actually
%                 executes the mesoRD command.
% ML 2015-09-10

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_runMesoRD.m, part of the SMeagol package
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
    error(['Cannot find MesoRD in the correct place : ' mesoRDstring])
end

if(~exist('targetFolder','var') || ~exist(targetFolder,'dir'))
    % then use the current folder
    targetFolder=pwd;
end
% hack to work around the use of pathnames with spaces: add " around the
% path to mesord:
mesoRDcommand=['"' mesoRDstring '" ' mesoRDoptions];
cd(targetFolder)

disp('---------------------')
disp(' full mesord command:')
disp(mesoRDcommand)
disp('---------------------')
try
    sysflag=system(mesoRDcommand);
catch me
    cd(startFolder)
    rethrow me
end
cd(startFolder)





