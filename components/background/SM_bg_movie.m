% IMBG=SM_bg_movie(t,optBG,cameraOpt)
% 
% A wrapper to use a movie as background, given as a tif stack. If the
% number of simulated frames exceed the number of frames in the background
% movie, the background movie starts over again. Since this is intended for
% use with actual experimental movies, no additional noise is added.
%
% input:
% t       : Start of frame for which background is to be generated  (not
%           used here). 
% optBG*  : background options struct (opt.background in Palantir).
%           Fields: 
% optBG.inputMovie     : path to input image file
%
% cameraOpt*: Palantir camera options struct (opt.camera), to pass
%             information about the image size and camera noise parameters.
%             
% * Only at first call, or when changing options. Otherwise, old values are
% remembered. When passing options, the internal counter is reset, and no
% movie is returned, so that the next call with time only returns the first
% background frame.
%
% output
% 1) IMBG   : with no input arguments, a default optBG struct is returned.
% 2) IMBG   : with at least the t input parameter, the next frame of the
% input movie is returned (independent of the value of t). When the movie
% runs out of frames, it is restarted.
% ML 2015-01-23

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_bg_movie, simulated image background in the SMeagol package
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
function [IMBG,frame]=SM_bg_movie(~,optBG,cameraOpt)

persistent bgTifStack stackCount stackSize;

if(nargin==0) % no input arguments, return default opt struct
    IMBG=struct('inputMovie',''); % default: empty string as path for input movie
    return
elseif(nargin==3) % then reinitialiation is called for
    stackCount=0;
    % read movie to use as background
    [PATHSTR,NAME,EXT] = fileparts(optBG.inputMovie);
    % SM_runsimulation uses double movie matrices, even though most tifs
    % are uint16.
    bgTifStack=ML_loadStack2(fullfile(PATHSTR,[NAME EXT]));
    clear PATHSTR NAME EXT

    % check size
    [ny,nx,stackSize]=size(bgTifStack);
    
    if(nx ~= cameraOpt.xrange_px || ny ~= cameraOpt.yrange_px)
        error('SM_bg_movie: movie size does not match camera ROI size.')
    end
    IMBG=[];
    frame=0;
    return
end
% update frame counter and extract the correct background frame
stackCount=mod(stackCount,stackSize)+1;
frame=stackCount;
IMBG=bgTifStack(:,:,frame);
