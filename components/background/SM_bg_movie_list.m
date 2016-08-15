% [IMBG,frame]=SM_bg_movie_list(~,optBG,cameraOpt)
% 
% A wrapper to use a list (cell vector of file names) of movie files as
% background, assuming that the formats match. If the number of simulated
% frames exceed the total number of simulated frames in the previous
% simulation, background movie, the background starts over. No additional
% random background is added, but offset and readout noise must be set to
% zero by hand in order to preserve the background faithfully. 
%
% input:
% t       : Start of frame for which background is to be generated  (not
%           used here). 
% optBG*  : background options struct (opt.background in Palantir).
%           Fields: 
% optBG.movieName : cell vector of path to input image files. To use a
%           previous SMeagol simulation, pass the movieName field of the
%           output structure from SM_getResults.
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
% ML 2016-02-20

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_bg_movie_list, simulated image background in the SMeagol package
% ========================================================================= 
% Copyright (C) 2016 Martin Lindén and Johan Elf
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
function [IMBG,frame]=SM_bg_movie_list(~,optBG,cameraOpt)

persistent bgTifStack frameCount stackSize movieName movieCount;

if(nargin==0) % no input arguments, return default opt struct
    IMBG=struct('movieName',{{'movie1.tif','movie2.tif','etc.tif'}}); % default: empty string as path for input movie
    return
elseif(nargin==3) % then reinitialiation is called for
    frameCount=0;
    movieName=optBG.movieName;
    movieCount=1;
    % read movie to use as background
    [PATHSTR,NAME,EXT] = fileparts(movieName{movieCount});
    % SM_runsimulation uses double movie matrices, even though most tifs
    % are uint16.
    bgTifStack=ML_loadStack2(fullfile(PATHSTR,[NAME EXT]));
    clear PATHSTR NAME EXT

    % check size
    [ny,nx,stackSize]=size(bgTifStack);
    
    if(nx ~= cameraOpt.xrange_px || ny ~= cameraOpt.yrange_px)
        error('SM_bg_movie_list: movie size does not match camera ROI size.')
    end
    IMBG=[];
    frame=0;
    return
end
% update frame counter and extract the correct background frame
frameCount=frameCount+1;
% switch to next movie if needed
if(frameCount>stackSize)
    frameCount=1; % restart frame counter
    movieCount=mod(movieCount,length(movieName))+1; % update movie counter 
    
    % read new movie
    [PATHSTR,NAME,EXT] = fileparts(movieName{movieCount});
    % SM_runsimulation uses double movie matrices, even though most tifs
    % are uint16.
    bgTifStack=ML_loadStack2(fullfile(PATHSTR,[NAME EXT]));
    disp(['SM_bg_movie_list loaded new movie ' fullfile(PATHSTR,[NAME EXT])])
    clear PATHSTR NAME EXT
end

IMBG=bgTifStack(:,:,frameCount);
