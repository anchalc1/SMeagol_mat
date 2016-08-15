% [trj,exitflag]=SM_parse_mesoRD_logfile(filename,lengthConversion,timeConversion)
%
% parse a mesoRD tracking logfile and generate SMeagol options in the form
% of the opt.trj substruct. 
%
% In addition to the standard fields, the units used in the MesoRD
% simulation is also reported, in the form of the fields
% D_input_unit
% voxelSize_input_unit
% voxelSize_input
% time_input_unit
%
% lengthConversion,timeConversion are conversion factors applied to the
% extracted values before rturning them, e.g., 
% lengthConversion=1e7, timeConversion=1/60
% converts cm -> nm, s -> min, cm^2/s -> nm^2/min.
% If not given, SM_parse_mesoRD_logfile will launch a small GUI to ask for
% them.
% 
% timeConversion is saved as trj.timeScale
%
% exitflag=0 after a successful conversion, and >0 otherwise. 
%
% Presently, up to 10000 species can be parsed. Edit the variable
% max_species in this file if more species are needed.
%
% ML 2015-03-02

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_license.m, prints a short license text for the SMeagol package
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
function [trj,exitflag]=SM_parse_mesoRD_logfile(filename,lengthConversion,timeConversion)

exitflag=1;
trj=struct;
graphicalmode=false;
max_species=10000;

if(nargin<3 || isempty(lengthConversion) || isempty(timeConversion) )
    % then additional input will be asked for, and we can go to
    % graphicalmode
    graphicalmode=true;
end 

% test for file existence
if(~exist(filename,'file'))
   if(graphicalmode)
       errordlg(['SM_parse_mesoRD_logfile could not find the mesoRD file ' filename])   
    return
   else
       error(['SM_parse_mesoRD_logfile could not find the mesoRD file ' filename])          
   end
end

% start parsing
fid=fopen(filename);

% discard list of species heading
tline=fgetl(fid);
% parse species list
trj.speciesNames={};

trj.D=[];
trj.D_input_unit='';

trj.degradedName='-1';  % name of the degraded state
trj.voxelSize_input=0; % length conversion factor from voxel to nm [nm/voxel]
%trj.voxelSize=0; % length conversion factor from voxel to nm [nm/voxel].
%This field is created later to make the output in the GUI easier to
%understand.
trj.voxelSize_input_unit='';
% Note that SMeagol no longer needs to know the timestep of the input
% trajectory; this is read off directly from the data instead.
trj.time_input_unit='';

for s=1:max_species
    tline=fgetl(fid);
    if(~ischar(tline)) % then end-of-file came early
        return
    end
    if(length(tline)>=22 && strcmp(tline(1:22),'Size of the subvolume:'))
        break % All species parsed!
    end
    ind=strfind(tline,' '); % typical line: FN 3e-08 cm2ps
    trj.speciesNames{s}=tline(1:ind(1)-1);
    trj.D(s)=str2double(tline(ind(1)+1:ind(2)-1));
    trj.D_input_unit=tline(ind(2)+1:end);
end
% voxel length
ind=strfind(tline,':')+1;
tline=tline(ind(1):end); % remove text before ':'
ind=strfind(tline,' ');
trj.voxelSize_input=str2double(tline(ind(1):ind(2)));
trj.voxelSize_input_unit=tline(ind(2)+1:end);

% time step
tline=fgetl(fid);
ind=strfind(tline,':')+1;
tline=tline(ind(1):end); % remove text before ':'
ind=strfind(tline,' ');
trj.time_input_unit=tline(ind(2)+1:end);

fclose(fid);

% unit conversion factors
if(graphicalmode)
    [lengthConversion,timeConversion]=mesoRD_unitconversionGUI(trj);
end
trj.timeScale=timeConversion;
trj.voxelSize=trj.voxelSize_input*lengthConversion;
trj.D=trj.D*lengthConversion^2/timeConversion;

exitflag=0;

% sample file
%List of species and their corresponding diffusion constant in their current compartment:
%FN 3e-08 cm2ps
%BN 4e-10 cm2ps
%FC 3e-08 cm2ps
%BC 0 cm2ps
%Size of the subvolume: 1e-06 cm
%Timestep between the outputs: 0.001 sec
