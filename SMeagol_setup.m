function SMeagol_setup()
% Add paths that are needed to run SMeagol to the matlab path.
%
% Also look for core/SM_mesoRD_tracking_path.mat where the path to the
% mesoRD binary is stored: if that file does not exist, run
% SM_find_mesoRD_path in the 'no userinput' mode, and issue a warning it
% the binary is not found.

dir0=pwd;
[SMeagolFolder,~,~]=fileparts(mfilename('fullpath'));

addpath(SMeagolFolder)
addpath(genpath([SMeagolFolder filesep '.' filesep 'components']))
addpath(genpath([SMeagolFolder filesep '.' filesep 'core']))
addpath(genpath([SMeagolFolder filesep '.' filesep 'gui']))
SM_license('simulated superresolution microscopy')
disp('----------------------------------------------------------------------')
disp('Added SMeagol paths. The main GUI is called SMeagol_gui.')
disp('----------------------------------------------------------------------')
cd(dir0);
