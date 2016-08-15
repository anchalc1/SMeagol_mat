% -------------------------------------------------------------------------
% SMeagol -- Simulated Single Molecule fluorescence Microscopy
% -------------------------------------------------------------------------
% Copyright (C) 2015 Martin Lindén and Johan Elf
%
% This program comes with ABSOLUTELY NO WARRANTY. This is free software,
% and you are welcome to redistribute it under certain conditions. See
% license.txt for details.
%
% Additional permission under GNU GPL version 3 section 7 
%
%   If you modify this Program, or any covered work, by linking or
%   combining it with Matlab or any Matlab toolbox, the licensors of this
%   Program grant you additional permission to convey the resulting work. 
%
% -------------------------------------------------------------------------
%
% High-level functions
% -------------------------------------------------------------------------
% SMeagol_setup   : Set up SMeagol paths etc. Add to startup.m to add 
%                   SMeagol to your matlab path at startup. 
% SMeagol_gui     : Start SMeagol GUI to edit parameters and run
%                   simulations.
% SM_runsimulation: run SMeagol simulation from command line.
% -------------------------------------------------------------------------
%
% Tests and configurations
% -------------------------------------------------------------------------
% SMeagol_setup     : Add SMeagol to the matlab path, and run a
%                     light-weight searhc for the mesoRD engine.
% SM_mesoRD_tracking_test : test that the integrated mesoRD works
% -------------------------------------------------------------------------
%
% Manipulate runinput files and options:
% -------------------------------------------------------------------------
% SMeagol_gui     : Start SMeagol GUI to edit parameters and run
%                   simulations. 
% SM_getOptions   : Read runinput file to options struct.
% SM_getResults   : Load the results from a finished simulation.
% SM_writeRuninputFile :  Write runinput file from options struct.
% SM_parse_mesoRD_logfile : Lowlevel tool to parse
%                           molecule_tracking_log.txt files from  MesoRD.
% SM_opt2str      : Lowlevel tool to convert options struct to string cell
%                   vector. 
% SM_str2opt      : Lowlevel tool to convert string cell vector to options
%                   struct. 
% -------------------------------------------------------------------------
%
% Tools for computing relative paths
% -------------------------------------------------------------------------
% SM_relative_path        : Computes the relative path between two
%                           absolute paths. 
% SM_relative_path_to_file: Computes the path to a file relative to a given
%                           starting path. 
% -------------------------------------------------------------------------
%
% Low-level simulation functions
% -------------------------------------------------------------------------
% SM_markovForward : Simulate a continuous time Markov process (used for
%                    the photophysics in SMeagol). 
% SM_xyToImage     : Convert positions of detected photons to an image.
% brownianBridge_piecewise : construct Brownian bridges.
