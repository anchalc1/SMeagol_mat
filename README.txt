 SMeagol is a program suite for simulating superresolution microscopy
 data with diffusing fluorophores, which should be useful for
 development and validation of methods in camera-based single particle
 tracking. The program runs on Matlab.
 https://github.com/bmelinden/SMeagol_mat is the matlab source code.
 
 The latest release, including install scripts, documentation and a
 forum for discussion and questions can be found at
 http://smeagol.sourceforge.net/
 =========================================================================
 Copyright (C) 2015 Martin Lindén and Johan Elf
 
 E-mail: bmelinden@gmail.com, johan.elf@gmail.com
 =========================================================================
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or any
 later version.  This program is distributed in the hope that it will
 be useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
 the GNU General Public License for more details.
 
 Additional permission under GNU GPL version 3 section 7
  
 If you modify this Program, or any covered work, by linking or
 combining it with Matlab or any Matlab toolbox, the licensors of this
 Program grant you additional permission to convey the resulting work.
 
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 =========================================================================
 If you use SMeagol for research, please cite our work: 
 Lindén, M., V. Ćurić, A. Boucharin, D. Fange, and J. Elf. Simulated single 
 molecule microscopy with SMeagol. Bioinformatics. 32: 2394-2395 (2016).

 SMeagol is written by Martin Lindén, David Fange, and Alexis
 Boucharin, and also includes software developed by Martin Lindén and
 Fredrik Persson. See copyright in individual files in external/.
=========================================================================
Release notes:
current:
 - changed behavior of SM_runsimulation.m to always generate
   opt.output.maxFrame frames even if there is not enough input
   trajectories. Exception: with opt.output.maxFrame=inf, the
   simulation is terminated when the input trajectories end.

v1.0.2 (2016-08-03): bugfix and small additions
 - new photoactivation model
 - new background models
 - bugfix: multi-movie spotStats was sorted incorrectly
 - bugfix: fixed typo that made SM_psf_GibsonLanni_680nm_NA149 actually 
   used  same PSF data as SM_psf_GibsonLanni_580nm_NA14
v1.0.1 (2015-12-11): bugfix and small additions
 - bugfix in mesoRD_unitconversionGUI
 - added decaying uniform background model
 - added single pulse activation model
v1.0  (2015-07-xx): beta release for review purposes. 
 
=========================================================================
Installation:

1. Automated full installation: Installation scripts including
pre-compiled binaries of the stochastic reaction-diffusion simulation
enging can be found on http://smeagol.sourceforge.net/

2. Matlab part and development: The Matlab part of SMeagol, available
at at https://github.com/bmelinden/SMeagol_mat, can run
independently. To add the reaction-diffusion part, one can simply copy
an appropriate mesord binary (from e.g., method 1) to
/path/to/SMeagol_mat/mesord/mesord

