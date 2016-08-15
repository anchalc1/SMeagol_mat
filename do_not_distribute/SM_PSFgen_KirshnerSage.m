% R=SM_PSFgen_KirshnerSage()
%
% Function to compute and parameterize a ROTATIONALLY SYMMETRIC 3D PSF
% model based on the PSF generator by Kirshner & Sage [1,2], available at
% http://bigwww.epfl.ch/algorithms/psfgenerator/  
% 
% Before running this function, download the matlab plugin of the
% PSFGenerator and put it in your matlabpath as described (or do it using
% this function).
%
% Standard usage with no input
% 
% 1) Add the PSFGenerator to the matlabpath if not already present.
% 2) The PSFGenerator GUI is started. Choose model and setting. Hints: 
%   - the computational focal plane will be placed in the middle of the
%     z-stack, which agrees with the simulated PSF is 'Particle Psoition Z'
%     is set to zero.
%   - Use high enough resolution compared to the pixels you will use for
%     simulation and the features of the PSF.
%   - this calculation only requires the PSF in the plane x=0. Hence,
%       a) use an odd number of pixels in the x-direction, and 
%       b) use the minimum number (five) of pixels in the x-direction to
%       save memory. 
%   - The y- and z-ranges should be large enough to capture most of PSF.
%     Photons outside the simulated range will be sent off to infinity.
% 3) Select 'compute PSF' to simulate the current PSF settings
% and read it into matlab. Update PSFGenerator settings and recompute as
% necessary. (On java memory problems, try restarting matlab...).
% 4) When a good PSF has been computed, select 'parameterize'. Matlab
% will now aslk for the Pixelsize XY and Z-step parameter values, and then
% generate data for a lookup table representation of the last computed PSF,
% using the x=0 plane and the assumption that the computed PSF is
% rotationally symmetric around the z-axis.
% 5) Select a .mat file to save the PSF parameters in. To simulate the PSF
% in Palantir, give this file as input parameter to the
% SM_psf_numericalPSF_Ruzp psf model.
%
% references (see http://bigwww.epfl.ch/algorithms/psfgenerator/#ref):
% 1.    Hagai Kirshner, François Aguet, Daniel Sage, Michael Unser, 3-D PSF
%       Fitting for Fluorescence Microscopy: Implementation and
%       Localization Application, Journal of Microscopy, vol. 249, no. 1,
%       pp. 13-25, 2013.
% 2.    Alessandra Griffa, Nathalie Garin and Daniel Sage, Comparison of
%       Deconvolution Software in 3D Microscopy. A User Point of View, Part
%       I and Part II, G.I.T. Imaging & Microscopy, vol 1, pp. 43-45, 2010.
%
% ML 2015-01-19

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_PSFgen_KirshnerSage, PSF-constructor in the Palantir package
% =========================================================================
% Copyright (C) 2015 Martin Lindén
%
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%
%  If you modify this Program, or any covered work, by linking or combining
%  it with Matlab or any Matlab toolbox, the licensors of this Program
%  grant you additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% actual code
function psfData=SM_PSFgen_KirshnerSage(varargin)
%% input parameters
if(nargin>0)
    warning('SM_PSFgen_KirshnerSage does not currently support input parameters.')
end
%% PSFGenerator.jar in the java path?
try
    PG=PSFGenerator; % initialize a PSFGenerator
catch me
    
    [PSFfile,PSFpath]=uigetfile({'PSFGenerator.jar','PSFGenerator matlab plugin'},...
        'Select the PSFGenerator matlab plugin .jar file.');
    if(strcmp(PSFfile,'PSFGenerator.jar'))
        javaaddpath(fullfile(PSFpath,PSFfile));    
        PG=PSFGenerator; % initialize a PSFGenerator
    else
       errordlg('PSFGenerator.jar matlab plugin not found. Quitting.')
       R=struct;
       return
    end
end
%% open PSF gui and wait for user confirmation of good parameters
PG.gui;
psf=[]; % psf matrix
% set up dialog box for the y and z lattice spacing
Xstep=100;
Zstep=250;
while(true)
    PSFaction=questdlg_normalWindow('Select PSF action:','PSF menu',...
        'compute PSF','parameterize','quit','compute PSF');
    switch PSFaction
        case 'quit'
            psfData=struct;
            return
        case 'compute PSF'
            % ask for some parameter values
            PSFpar=inputdlg({'Pixelsize XY:','Z-step:'},...
                'Please enter PSF parameters',1,{num2str(Xstep),num2str(Zstep)});
            if(isempty(PSFpar))
                continue
            end
            Xstep=str2double(PSFpar(1));
            Zstep=str2double(PSFpar(2));
            if (isnan(Xstep) || isnan(Zstep) )
                errordlg('Unreadable PSF parameters. Not computing.')
                continue
            end
            % compute PSF
            psf=PG.get; 
            % visualize the XZ plane and radial distribution
            % note interchanged xy-indices in the psf matrix
            pltPSF_yz_rz(psf,Xstep,Zstep,11);
            
        case 'parameterize'
            % make some basic tests, and go on to parameterization if pass
            if(isempty(psf))
                errordlg('Empty PSF. Compute it first.')
                continue;
            elseif(size(psf,2)/2 == round(size(psf,2)/2) ) 
                % check for odd x-size (so that x=0 is present)
                % note interchanged xy-indices in the psf matrix
                errordlg('Odd number of X pixels needed. Set first entry of Size XYZ to a (low) odd number and recompute.')
                continue
            else
                break
            end
    end
end
%% parameterize
[nx,ny,nz]=size(psf); % psf size: note interchanged xy-convention ued by PSFGenerator
x_range=Xstep*((1:nx)-mean(1:nx));
r_range=x_range(x_range>=0);
z_range=Zstep*(((1:nz)-mean(1:nz)));

% cumulative normalization for every z
Fcum=zeros(ceil(nx/2),nz);
PfiniteR=zeros(1,nz);

parfor kz=1:nz
    tic
    psfY1Z=reshape(psf(:,ceil(ny/2),kz),1,nx); % this particular slice (y=0,z fixed)
    
    Px=spline(x_range,psfY1Z); % spline interpolation of psf(x=0,y,z=z_k);
    
    dFdr=@(t,f)(2*pi*t.*ppval(Px,t).*ones(size(f)));
    opt=odeset('reltol',1e-9,'abstol',1e-9);
    [~,F]=ode45(dFdr,r_range,0,opt);
    PfiniteR(kz)=F(end);     
    Fcum(:,kz)=F/F(end); % conditional cumulative radial distribution         
    disp([int2str([kz nz]) ' ' num2str(toc)])    
end
PfiniteR=PfiniteR/max(PfiniteR); % probability of finite radius
% invert the conditional cumulative distribution
u=linspace(0,1,500);        % U(0.1) random number
R_uz=zeros(length(u),nz);    % corresponding radius, conditional on it being finite
parfor kz=1:nz
    R_uz(:,kz)=interp1(Fcum(:,kz),r_range,u,'spline');
end

% test constructing a gridded interpolant
FRuz=griddedInterpolant({u,z_range},R_uz,'spline');
%% save interpolation raw data
% construct data struct
psfData=struct('R_uz',R_uz,'u',u,'z',z_range,'PfiniteR',PfiniteR);

% choose save file
while(true)
    [PSFfile,PSFpath]=uiputfile({'*.mat','matlab .mat file'},'Save PSF parameters.');
    if(isnumeric(PSFpath) && isnumeric(PSFfile))
        % then no savefile was selected
        PSFaction=questdlg_normalWindow('No output file selected. Quit without saving?',...
            '','yes','no','no');
        switch PSFaction
            case 'yes'
                return
            case 'no'
                continue
        end
    end
    break
end
save(fullfile(PSFpath,PSFfile),'-struct','psfData');
end
function pltPSF_yz_rz(psfMatrix,dX,dZ,figureNumber)
%% plot the psf in various ways
[nx,ny,nz]=size(psfMatrix);
%psfYZ=reshape(psfMatrix(ceil(nx/2),:,:),ny,nz); % center slice
psfXZ=reshape(psfMatrix(:,ceil(ny/2),:,:),nx,nz); % center slice, unintuitive xy convention of PSFGenerator

R=dX*(((1:nx)-mean(1:nx))'*ones(1,nz));
X=dX*(((1:nx)-mean(1:nx))'*ones(1,nz));
Z=dZ*ones(nx,1)*(((1:nz)-mean(1:nz)));

figure(figureNumber)
clf

subplot(1,2,1)
plot(0,0,'k.')
surf(X,Z,psfXZ,'edgecol','none')
view([0 90])
xlabel('y')
ylabel('z')
axis equal
axis tight
title('PSF density in yz plane')


subplot(1,2,2)
plot(0,0,'k.')
surf(R,Z,psfXZ.*abs(R),'edgecol','none')
view([0 90])
axis equal
axis tight

xlabel('r')
ylabel('z')
title('PSF radial density')

end
