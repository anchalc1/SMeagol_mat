% vxl=SMRD_mesoRDgeometry(filename)
%
% read voxel midpoints from a mesoRD geometry file, where each row 
% contains the voxel coordinates and then a compartment name string.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMRD_mesoRDgeometry.m, part of the SMeagol package
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
function vxl=SMRD_mesoRDgeometry(filename)
% read mesoRD geometry file
% each line contains three coordinates and then a compartment name
fid=fopen(filename,'r');
numLines=0;
disp(['SMRD_mesoRDgeometry parses ' filename ' ...'])
vxl=zeros(3,1e5);
tic
while(true)
    linestr=fgetl(fid);
    if ~ischar(linestr), break, end
    numLines=numLines+1;
    % identify the first three numbers
    spcind=strfind(linestr,' '); % find all spaces on the line
    spcind=spcind(1+find(diff(spcind)>1)); % space after each word
    vxl(1:3,numLines)=str2num(linestr(1:spcind(3)));
end
fclose(fid);
vxl=vxl(1:3,1:numLines);
disp(['SMRD_mesoRDgeometry finished ' int2str(numLines) ' voxels in '  filename ' in ' num2str(toc) ' s.'])
clear spcind linestr fid numLines;

