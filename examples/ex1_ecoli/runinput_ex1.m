% SMeagol runinputfile, created 08-Sep-2015 13:59:06.
% a note on units: SMeagol does not know about units, and so the user is
% charged with using consistent units of length, time, and diffusion
% constants (units of length^2/time).
%% trj
% ----------------------------------------------------------------------- %
% trj: information about the input data.
trj.reactionFile='./mesoRD_2state_noturnover/reactions.txt';
trj.trajectoryFile='./mesoRD_2state_noturnover/trajectories.txt';
trj.degradedName='-1';
trj.timeScale=1;
trj.voxelSize=10;
trj.speciesNames{1}='FN';
trj.speciesNames{2}='BN';
trj.speciesNames{3}='FC';
trj.speciesNames{4}='BC';
trj.D=[ 1e+07       1e+05       1e+07           0];
% ----------------------------------------------------------------------- %
%% output
% ----------------------------------------------------------------------- %
% output: what to output, and where.
output.resultFile='./results/ex1_ecoli_results.mat';
output.writeTifMovie=true;
output.plotTifMovie=true;
output.showPhotons=false;
output.showEmitters=false;
output.plotTrj=false;
output.movieLength=100;
output.maxFrames=200;
output.plotTrjZRange=[ -500         500];
output.movieFormat='tiff';
output.movieOptsImwrite = {};
% ----------------------------------------------------------------------- %
%% sample
% ----------------------------------------------------------------------- %
% sample: parameters describing illumination and image capture. SMeagol
% basically assumes that illumination and aquisition coincide, but
% continuous illumination can be modeled by an appropriate choice of
% photophysics-parameters.
sample.dt=0.009;
sample.tE=0.001;
% ----------------------------------------------------------------------- %
%% activation
% ----------------------------------------------------------------------- %
% activation: parameters describing the fluorophore activation process.
activation.t1=0;
activation.ta=0.1;
activation.Pa=0.3;
activation.td=0;
activation.ka=[ 0           0];
activation.type='SM_photoActivation_Pa_basal';
% ----------------------------------------------------------------------- %
%% baseIntensity
% ----------------------------------------------------------------------- %
% baseIntensity: every fluorescent group in the simulation has a basic
% emission intensity (photons/time) during illumination, which can vary
% from molecule to molecule, as determined at activation by these
% parameters.
baseIntensity.intensity=2e+05;
baseIntensity.type='SM_activationIntensity_uniform';
% ----------------------------------------------------------------------- %
%% photophys
% ----------------------------------------------------------------------- %
% photophys: parameters describing the dynamics of blinking and bleaching
% in terms of a Markov process (independent of the diffusive states
% described by the input trajectories).
photophys.bleach_time=0.5;
photophys.type='SM_fluo_only_bleach';
% ----------------------------------------------------------------------- %
%% psf
% ----------------------------------------------------------------------- %
% psf: parameters to simulate the microscope point-spread-function, i.e.,
% the (stochastic) map from the position of a fluorophore as it emits a
% photon to the position on the camera chip at which that photon is
% detected.
psf.type='SM_psf_GibsonLanni_584nm_NA14';
% ----------------------------------------------------------------------- %
%% camera
% ----------------------------------------------------------------------- %
% camera: these parameters describe a) The region of interest (ROI), i.e.,
% the size, shape, and location of the region imaged by the camera, and b)
% the noise properties of the EMCCD chip.
% (a) is described in terms of the size and number of active camera pixels,
% plus a linear transformation of simulated coordinates x to
% camera-centered coordinates y, given by y = (voxelsize)*A*x+b, where A is
% a 3*3 matrix, and b is a 3*1 vector. 
% (b) is parameterized in terms of the offset, readout noise (standard
% deviation), and EM gain (average number of photons per camera count). We
% use the model described in the Mortensen et al. (Nat Meth 7, 377–381,
% 2010, doi: 10.1038/nmeth.1447).
camera.alpha=0.025;
camera.offset=100;
camera.sigmaReadout=4;
camera.pixLength=80;
camera.xrange_px=35;
camera.yrange_px=65;
camera.A(1,:)=[0.22495    -0.97437           0];
camera.A(2,:)=[0.97437     0.22495           0];
camera.A(3,:)=[0           0           1];
camera.b=[ 1471.6      2700.7        1.19]'; %(note transpose!)
% ----------------------------------------------------------------------- %
%% background
% ----------------------------------------------------------------------- %
% background: parameters to describe how the noisy image background is to
% be generated.
background.photons_per_pixel=1;
background.type='SM_bg_constant';
% ----------------------------------------------------------------------- %
