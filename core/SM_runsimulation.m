% SM_runsimulation : run SMeagol microscopy simulations
%
% function opt=SM_runsimulation(runinput)
%
% Runs a SMeagol micrsocopy simulation using the parameters specified in
% runinput.
%
% runinput  : The name of the runinput file, or a SMeagol options
% struct. If the latter, SM_runsimulation also looks for the field
% rStateInit and uses it to reset the randstram state with the rng command.
%
% output:
% opt   - an options struct with the parameters used in the simulation,
%         plus two extra fields. 
%           - rStateInit: the randstream state. If given in the input opt
%             structure, the randstream will be reset, and thus it is
%             possible to rerun a simulation with the same random numbers.  
%           - runVariables.runID : a string with date and time among other
%             things. Useful to keep track of when something was run.
%         The output struct is also saved in the file specified in the
%         output.resultFile parameter.
%
% Results file: the file specified in the output.resultFile parameter also
% contains a few fields that can be used as 'ground truth' to test various
% analysis methods.
%   - positionSnapshots: fluorescent particle positions at the beginning of
%     each frame.
%   - emissionAverage  : average position of photon emission events for
%     each fluorescent particle during each frame. 
%   - detectionAverage : average position of detected photons from every
%     fluorescent particle in each frame. 
%   - diffusiveState: hidden diffusive state of each fluorescent particle
%     in the beginning of each frame.
%   - spotStats : some information about the emitted photons for each spot.
%   - movieName : a map from movieNr (below) to tif files (if saved).
% Layout of the ground truth fields: one cell element for each particle.
% The rows are as follows: particle k, located in frame frameNr in
% movieName{movieNr}.
%   positionSnapshots{k}= [x(t) y(t) z(t) frameNr movieNr]; 
%   emissionAverage{k}  = [x(t) y(t) z(t) frameNr movieNr];
%   detectionAverage{k} = [x(t) y(t) 0    frameNr movieNr];
%       (x(t),y(t) = NaN if no photons where emitted).
%   diffusiveState{k}   = [s(t) frameNr movieNr]; s = hidden state
%   spotStats{k}        = [nPhotons covD(1,1) covD(1,2) covD(2,2)];
%   covD is the covariance matric of the detected photon positions.
%
% ML 2015-04-28


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_runsimulation.m, main simulation engine of the SMeagol package
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

function opt=SM_runsimulation(runinput)

SM_license('SM_runsimulation')
%% read and check options

if(isstruct(runinput)) % then opt struct is passed. Rerun with same random seed if possible
    opt=runinput;
    %runiput=fullfile(opt.runinputroot,opt.runinputfile); % just to do something with it
    
    if(isfield(opt,'rStateInit'))
        rng(opt.rStateInit)
    else
       disp('No random seed opt.rStateInit found in input struct. Exact reproduction not possible.')
    end
else % runinput is a runinput file name
    opt=SM_getOptions(runinput);
    opt.rStateInit=rng; % save random generator state for reproducibility
end
% generate a run-id
opt.runVariables.runID=[opt.runinputfile ' : ' datestr(now)];

% save random seed to results file
% create results file name & folder
[mpath,mname,~]=fileparts(opt.output.resultFile);
fullResultName=fullfile(opt.runinputroot,mpath,mname);
[mk_flag,mk_message]=mkdir(opt.runinputroot,mpath);
disp(['runinput folder : ' opt.runinputroot])
disp(['result file     : ' opt.output.resultFile])
if(~mk_flag && ~isempty(mpath))
    disp('Could not create target folder for results file.')
    error(mk_message);
else
    disp('saved simulation options to result file.  Starting simulation.')
    disp('--------------------------------------------------------------')
end
save(fullResultName,'opt');

% exposure time must be smaller than sampling time
if(opt.sample.tE>opt.sample.dt)
    error('SMeagol: sample.tE > sample.dt not allowed. Exposure time must be smaller than sampling time.')
elseif(opt.sample.tE==opt.sample.dt)
    disp( 'SMeagol: sample.tE = sample.dt. Adjusting exposure time slightly lower to make tE < dt.')
    opt.sample.tE=opt.sample.tE*(1-10*eps);
end
% add degraded state to species names
allSpecies=[opt.trj.speciesNames {opt.trj.degradedName}];
%% read trajectory information
fiT=fopen(fullfile(opt.runinputroot,opt.trj.trajectoryFile),'r');
if(fiT==-1)
    error(['could not open trajectory file ' opt.trj.trajectoryFile ])
end
fiR=fopen(fullfile(opt.runinputroot,opt.trj.reactionFile),'r');
if(fiR==-1)
    error(['could not open reaction file ' opt.trj.trajectoryFile ])
end
%% initialize
% activation
actFunHandle=eval(['@' opt.activation.type]);
actFun=@(x)(actFunHandle(x,opt.activation,opt.sample));

% base intensity at activation 
baseIntFunHandle=eval(['@' opt.baseIntensity.type]);
baseIntFun=@()(baseIntFunHandle(opt.baseIntensity));

% pointspread function
psfFunHandle=eval(['@' opt.psf.type]);
psfFun=@(t,x)(psfFunHandle(t,x,opt.psf));

% photophysics (blinking and bleaching)
% pre-compute photophysics parameters
photophys=feval(opt.photophys.type,opt.photophys); 

% image background
bgFunHandle=eval(['@' opt.background.type]);
bgFunHandle(0,opt.background,opt.camera); % initialize
bgFun=@(t)(bgFunHandle(t));               % subsequent calls only need to supply frame time

% initialize SM_xyToImage
SM_xyToImage([0 0],opt.camera);

% movie number for multiple tif outputs
movieNumber=0; % gets updated to 1 in the first pasre-iteration
result_movieName=cell(0); % cell vectopr of move file names

% initialize plots
tifAxis=[];
trjAxis=[];
if(opt.output.plotTifMovie && opt.output.plotTrj) % display movie and trj in two subplots    
    figure(1)
    clf    
    tifAxis=subplot(2,1,1);
    trjAxis=subplot(2,1,2);
elseif(opt.output.plotTifMovie && ~opt.output.plotTrj) % display movie only 
    figure(1)
    clf    
    tifAxis=gca;
elseif(~opt.output.plotTifMovie && opt.output.plotTrj) % display trajectory only 
    figure(1)
    clf    
    trjAxis=gca;
end
% conversion of positions to camera frame
camA=opt.trj.voxelSize*opt.camera.A;
camB=opt.camera.b;
% conversion of input time
timeScale=opt.trj.timeScale;

% figure out the time step 
fiT=fopen(fullfile(opt.runinputroot,opt.trj.trajectoryFile),'r');
[t1,~,~,~]=SM_readLine(fiT,camA,camB,timeScale,allSpecies); % read initial configuration
[t2,~,~,~]=SM_readLine(fiT,camA,camB,timeScale,allSpecies); % read initial configuration
timestep=t2-t1;
fclose(fiT);

% set up initial configuration
fiT=fopen(fullfile(opt.runinputroot,opt.trj.trajectoryFile),'r');
[tTrj,idTrj,sTrj,xTrj]=SM_readLine(fiT,camA,camB,timeScale,allSpecies); % read initial configuration

% time margin for ignoring particles in the book-keeping
ignoreHorizon=3*max(timestep,opt.sample.dt); 

% initialize activation and trajectories
traj=cell(0);
activationTSI=zeros(3,max(idTrj));
%ignoreThis=zeros(1,size(activationTSI,2));
for kk=1:length(idTrj)
    k=idTrj(kk);
    traj{k}=[tTrj sTrj(kk) xTrj(kk,:)];
    activationTSI(1,k)=actFun(tTrj);
    activationTSI(2,k)=1;
    activationTSI(3,k)=baseIntFun();% basal fluorescence intensity of this fluorophore
end
clear k

tStart=tTrj;
nFrame=0;
%tFrame0=tStart+(nFrame-1)*opt.sample.dt; % start of this frame, redefined before first use
tFrame1=tStart+nFrame*opt.sample.dt; % end of this frame, initialization

tRct=tTrj;

% data structure for tracking truth
tru=struct;
tru.positionSnapshots=cell(0);
tru.emissionAverage=cell(0);
tru.detectionAverage=cell(0);
tru.diffusiveState=cell(0);
tic
%% parse
while(nFrame<opt.output.maxFrames)    
    %% reset frame variables
    nFrame=nFrame+1; % nFrame : 1 -> opt.output.maxFrames
    tFrame0=tFrame1;                     % start of this frame
    tFrame1=tStart+nFrame*opt.sample.dt; % end of this frame
    %%%traj0=traj;act0=activationTSI;tTrj0=tTrj;idTRJ0=idTrj;sTrj0=sTrj;xTrj0=xTrj;rState0=rng;      %%% debug variables
    % update movie number
    if(mod(nFrame,opt.output.movieLength)==1 || opt.output.movieLength==1)
        movieNumber=movieNumber+1;
    end            
    %% advance trajectory data past current frame
    ignoreThis=((activationTSI(2,:)==-1) ...                             % ignore degraded and bleached fluorophores
        +(activationTSI(1,:) > (tFrame0+ignoreHorizon)))>0;% and fluorophores that will not activate soon
    while(tTrj<=tFrame1)
        [tTrj,idTrj,sTrj,xTrj]=SM_readLine(fiT,camA,camB,timeScale,allSpecies,false,ignoreThis);
        if(isempty(tTrj))
            disp('no new line to be read (1)')
            break
        end
        
        for k=1:length(idTrj)
            if( idTrj(k)>length(traj))
                % new particle detected!
                traj{idTrj(k)}=[tTrj sTrj(k) xTrj(k,:)];
                activationTSI(1,idTrj(k))=actFun(tTrj);
                activationTSI(2,idTrj(k))=1;                
                activationTSI(3,idTrj(k))=baseIntFun();% basal fluorescence intensity of this fluorophore
            else
                traj{idTrj(k)}(end+1,:)=[tTrj sTrj(k) xTrj(k,:)];
            end
        end
        clear k
    end
    if(isempty(tTrj))
        disp('no new line to be read (2)')
        break
    end
    %% advance reactions past trajectory time tTrj
    % traj1=traj; %%% for debugging
    %act1=activationTSI;
    %tTrj1=tTrj;
    %idTRJ1=idTrj;
    %sTrj1=sTrj;
    %xTrj1=xTrj;
    %rState1=rng;  %%% for debugging
    while(tRct<=tTrj)
        [tRct,idRct,sRct,xRct]=SM_readLine(fiR,camA,camB,timeScale,allSpecies,true); % read all particles
        if(isempty(tRct)) % then there are no more rections to read
            break
        end
        % initialize if new particle
        if( idRct>length(traj) || activationTSI(2,idRct)==0 )
            % then this is a new particle, and needs an activation time
            activationTSI(1,idRct)=actFun(tRct);
            activationTSI(2,idRct)=1;
            activationTSI(3,idRct)=baseIntFun();% basal fluorescence intensity of this fluorophore
            % see if position needs to be saved
            if( activationTSI(1,idRct) > (tRct+ignoreHorizon) )
                ignoreThis(idRct)=true;
                traj{idRct}=[];
            else
                ignoreThis(idRct)=false;
                traj{idRct}=[tRct sRct xRct];
            end
        elseif(activationTSI(2,idRct)==-1) % this particle is already bleached or degraded, ignore it
            
        elseif(sRct==-1 && tRct < activationTSI(1,idRct)) % then degradation happened before activation
            activationTSI(1:2,idRct)=[tRct;-1];
            traj{idRct}=[];
            ignoreThis(idRct)=1;
        elseif( sRct > -1 && (tRct+ignoreHorizon) < activationTSI(1,idRct) )
            % activation is too far off at reaction time to take note.
            % coordinate will be read before activation anyway
        else
            traj{idRct}(end+1,:)=[tRct sRct xRct];
        end
    end
    %% sort events for each particle, and interpolate positions for t=tFrame1
    % traj2=traj;act2=activationTSI;rState2=rng; %%% for debugging
    for k=1:length(traj)
        if(~isempty(traj{k}))
            [~,ind]=sort(traj{k}(:,1));
            traj{k}=traj{k}(ind,:);
            
            ind=find(traj{k}(:,1)==tFrame1,1);
            if(isempty(ind)) % then interpolation to tFrame1 is needed
                i1=find(traj{k}(:,1)>tFrame1,1); % i1 is the first entry after tFrame1
                if(~isempty(i1) && i1>1)         % then interpolation is possible
                    i0=i1-1;                     % i0 is the last entry before tFrame1
                    s0=traj{k}(i0,2);
                    T=traj{k}(i0:i1,1);
                    X=traj{k}(i0:i1,3:5);
                    D=opt.trj.D(s0);
                    XI=brownianBridge_piecewise(T,X,D,tFrame1);
                    traj{k}=[traj{k}(1:i0,:);
                             tFrame1 s0 XI;
                             traj{k}(i1:end,:)];
                else % interpolation not possible
                    % In this case, the trajectory is simply not appended,
                    % and the particle was either just created or is about
                    % to degrade.
                    % disp('could not interpolate') %%% debug check
                    % keyboard
                end
            end
        end
    end
    clear s0 T X D XI k
    %% prune information for t<tFrame0-opt.sample.dt     % note that if the simulated time points are more sparse than this,    % points will be interpolated into the traj{} records anyway, so this    % pruning is safe irrespective of the spacing of simulated positions.
    % traj3=traj;act3=activationTSI;rState3=rng; %%% for debugging
    for k=1:length(traj)
        if(~isempty(traj{k}))
            % remove old events
            % ind=find(traj{k}(:,1)>=tFrame0-opt.sample.dt);traj{k}=traj{k}(ind,:);
            traj{k}=traj{k}(traj{k}(:,1)>=tFrame0-opt.sample.dt,:);
            % check for degradation clashes
            ind=find(traj{k}(:,2)==-1);
            if(~isempty(ind))
                tDegrade=traj{k}(ind,1);
                traj{k}=[traj{k}(traj{k}(:,1)<tDegrade,:);
                         traj{k}(ind,:)]; % keep only degradation, and the events strictly before it                
            end
        end
    end
    clear k
    %% simulate photophysics
    %%%traj4=traj;    act4=activationTSI;    rState4=rng;  %%% for debugging only
    xDetected=[]; % positions of detected photons in this frame, for plotting purposes
    xEmitted=[];  % positions of photon emission events in this frame, for plotting purposes
    for k=1:size(activationTSI,2)
        if(activationTSI(1,k)<tFrame1 && activationTSI(2,k)>0) % then this fluorophore is active in this frame
            try
                % states: 1 = fluorescent, 0=bleached
                t0=activationTSI(1,k); % current start time
                s0=activationTSI(2,k); % current photophysics state

                % advance through exposure time and emit photons
                if(t0<tFrame0+opt.sample.tE)
                    Q1 =photophys.Q{1};
                    QB1=photophys.kb{1};
                    
                    [~,sH]=SM_markovForward(Q1,s0,t0,tFrame0+opt.sample.tE,QB1);
                    % handle fluorophore degradation during the exposure time
                    if( traj{k}(end,2)==-1 && traj{k}(end,1)<tFrame0+opt.sample.tE)
                        %ind=find(sH(:,1)<traj{k}(end,1));  % keep only photoevents before degradation
                        %sH=[sH(ind,:) ;traj{k}(end,1) -1]; % terminate flourescence at degradation time
                        sH=[sH(sH(:,1)<traj{k}(end,1),:) ;traj{k}(end,1) -1]; % terminate flourescence at degradation time
                    end
                    activationTSI(1:2,k)=sH(end,1:2)';
                    %startingTrueState=s0; % remember true diffusive state for ground truth
                    startingTrueState=traj{k}(find(traj{k}(:,1)<=t0,1,'last'),2); % true diffusive state at the onset of exposure
                    %% now emit some light
                    ind=find(sH(1:end-1,2)>0); % starts of non-bleached periods
                    emF=reshape(photophys.emissionFactor(sH(ind,2)),length(ind),1);
                    tOn =sH(ind,1);
                    tOff=sH(ind+1,1);
                    dt = tOff-tOn; % fluorescent intervals
                    nMean=activationTSI(3,k)*dt.*emF; % average number of photons
                    nPh=poissrnd(nMean);         % emitted photons
                    
                    tPhk=[]; % photon emission times
                    for mm=1:length(nPh)
                        tPhk=[tPhk; tOn(mm)+dt(mm)*sort(rand(nPh(mm),1))];
                    end
                    if(~isempty(tPhk))
                        %tPhk=sort(tPhk);
                        % photon emission positions
                        T=traj{k}(:,1); %% error on line 243, index exceeds matrix dimensions?
                        D=opt.trj.D(traj{k}(1:end-1,2)); % last diffusion constant is not used anyway
                        X=traj{k}(:,3:5);
                        
                        % interpolate positions using Brownian bridges
                        xPhk=brownianBridge_piecewise(T,X,D,tPhk); % photon emission positions
                        yPhk=psfFun(tPhk,xPhk);                    % detected photon positions
                        
                        % detected photon positions
                        xEmitted=[xEmitted;   xPhk];
                        xDetected=[xDetected; yPhk];
                        
                        % compute and save tracking truth for this fluorophore
                        
                        %%% this is the position trace bug: traj{k}(1,:) is
                        %%% not always the beginning of the frame. Need to
                        %%% identify explicitly which row of traj{k} to use
                        %%% for x0.
                        [~,ix0]=min(abs(T-tFrame0));
                        %ix0=find(T>=tFrame0,1) % about twice as fast, but perhaps less robust
                        x0 =traj{k}(ix0,3:5); % position at beginning of frame
                        
                        
                        
                        xE =mean(xPhk,1);  % avergare photon emission position
                        xD =mean(yPhk,1);  % average position of detected photons
                        covD=cov(yPhk,1);  % covariance matrix of detected photons
                        if(numel(covD)~=4) % then probably no photons were emitted
                           covD=zeros(2,2); 
                           xD=[NaN NaN];
                        end
                        
                        % FP tracking format (we ignore cellNr and datasetNr).
                        % trajs{i} = [x(t) y(t) |xy(t)-xy(t-1)| frameNr cellNr datasetNr];
                        % we will save z(t) as well, but no step lengths:
                        % truTraj{i}      =[x(t) y(t) z(t) frameNr movieNr] %
                        % em/detAverage{i}=[x(t) y(t) z(t) frameNr movieNr] 
                        % z(t)=0 for the detectionAverage
                        % spotStats{i}=[nPh covXX covXY covYY];
                        if(length(tru.positionSnapshots)<k || isempty(tru.positionSnapshots{k})) % first frame in a trajectory
                            tru.positionSnapshots{k}         =[x0   nFrame movieNumber];
                            tru.emissionAverage{k}           =[xE   nFrame movieNumber];
                            tru.diffusiveState{k}            =[startingTrueState nFrame movieNumber];
                            tru.detectionAverage{k}          =[xD 0 nFrame movieNumber];
                            tru.spotStats{k}                 =[sum(nPh) covD(1,1) covD(1,2) covD(2,2)];
                        else
                            tru.positionSnapshots{k}(end+1,:)=[x0   nFrame movieNumber]; % add to existing trajectory
                            tru.emissionAverage{k}(end+1,:)  =[xE   nFrame movieNumber];
                            tru.diffusiveState{k}(end+1,:)   =[startingTrueState nFrame movieNumber];                            
                            tru.detectionAverage{k}(end+1,:) =[xD 0 nFrame movieNumber];
                            tru.spotStats{k}(end+1,:)        =[size(yPhk,1) covD(1,1) covD(1,2) covD(2,2)];
                        end
                    end
                    clear ind tOn tOff dt nMean nPh
                    % prepare for dark period
                end
                % advance until beginning of next frame
                s0=activationTSI(2,k); % current state
                t0=activationTSI(1,k); % current start time
                Q2=photophys.Q{2};
                QB2=photophys.kb{2};
                if(s0>0) % more Markov dynamics if the fluorophore is not bleached
                    [~,sH]=SM_markovForward(Q2,s0,t0,tFrame1,QB2);
                    activationTSI(1:2,k)=sH(end,1:2)';
                end
                % handle degradation during dark period (unless spontaneous
                % bleaching happened too).
                if( activationTSI(2,k)>-1 && traj{k}(end,2)==-1 && traj{k}(end,1)<tFrame1)
                    activationTSI(1:2,k)=[traj{k}(end,1); -1];
                end
            catch me
                disp('SM_runsimulation encountered an error during photophysics simulation')
                disp(me.message)
                disp(me.stack)
                errFile=['errorlog_SM_runsimulation' num2str(rand) '.mat'];
                save(errFile);
                disp(['saving workspace at crash to ' errFile ])
                rethrow(me)
            end            
        end
    end
    clear k
    %% make background
     BGcounts=bgFun(tFrame0);
    %% make image (if needed)
     if(opt.output.writeTifMovie || opt.output.plotTifMovie)
         IMphotons=SM_xyToImage(xDetected);         % detected photons from fluorophores
         IMcounts=round(gamrnd(IMphotons,1/opt.camera.alpha)); % EMCCD gain effect on detected photons
         IMemccd=BGcounts+IMcounts; % add signal and background
     end
    %% write to tif file and keep track of movie names 
    if(opt.output.writeTifMovie) 
        % create output filename with movie number in it
        [mpath,mname,~]=fileparts(opt.output.resultFile);
        % figure out the maximum number of digits in the movie number
        if(isfinite(opt.output.maxFrames))
            movieDigits=ceil(log10(opt.output.maxFrames/opt.output.movieLength+1));
        else
            movieDigits=8;
        end                
        mname_num_ext=sprintf(['%s_%0' int2str(movieDigits) 'd.tif'],mname,movieNumber);
        fullMovieName=fullfile(opt.runinputroot,mpath,mname_num_ext);
        result_movieName{movieNumber}=fullfile(mpath,mname_num_ext); % save path to this movie file
        
        [mk_flag,mk_message]=mkdir(opt.runinputroot,mpath);
            if(~mk_flag && ~isempty(mpath))
                disp(['runinput folder : ' opt.runinputroot])
                disp(['movie file      : ' result_movieName{movieNumber}])                
                disp('Could not create target folder for movie file.')
                error(mk_message);
            end
        
        % check for saturation
        Imax=max(max(IMemccd));
        if(Imax>65535)
            warning('image saturated, consider lowering the gain')
        end
        IMtif=im2uint16(uint16(IMemccd)); % convert to proper 16 bit image format
        
        if(mod(nFrame-1,opt.output.movieLength)==0)       % then this is the first frame in a new movie
            imwrite(IMtif,fullMovieName,...
                opt.output.movieFormat,'WriteMode','overwrite',opt.output.movieOptsImwrite{:});
        else
            imwrite(IMtif,fullMovieName,...
                opt.output.movieFormat,'WriteMode','append',opt.output.movieOptsImwrite{:});
        end
        if(mod(nFrame,100)==1 || nFrame==opt.output.maxFrames)
            disp(['wrote frame ' int2str(nFrame) ' to ' mname_num_ext ...
                ', ' num2str(toc/nFrame) ' s/frame.'])
        end
    end
    %% plot tif movie
    if(opt.output.plotTifMovie)
        if(~ishandle(tifAxis)) % then the axis was probably deleted
            errordlg('SM_runsimulation: movie window not found (closed?). Terminating simulation.')
            break
        end        
        
        axes(tifAxis)
        cla
        hold off
        
        % plot as an image and add detected photon positions on top to
        % match coordinate systems, use normal y-direction in order to
        % be comparable with the palantirAB_gui.
        imagesc(IMemccd);
        axis equal
        axis([0 opt.camera.xrange_px 0 opt.camera.yrange_px]+0.5)
        set(tifAxis,'Ydir','normal')
        %axis tight
        hold on
        if(opt.output.showPhotons && ~isempty(xDetected))
            plot((xDetected(:,1)/opt.camera.pixLength),...
                (xDetected(:,2)/opt.camera.pixLength),'r.','markersize',1)
        end
        if(opt.output.showEmitters && ~isempty(xEmitted))
            plot((xEmitted(:,1)/opt.camera.pixLength),...
                (xEmitted(:,2)/opt.camera.pixLength),'bx','markersize',1)
        end
        %grid on
        xlabel('x [px]')
        ylabel('y [px]')
        
        colormap gray
        title(['time: [ ' num2str([tFrame0 tFrame0+opt.sample.tE tFrame1]) ' ]'])
        colormap gray
        %set(gca,'clim',opt.output.tifCRange) % easier to just use imagesc
        %and have one less parameter
        
        % a dirty hack to write e.g., avi movies directly. Uncomment, and
        % create writerObj as described in doc videoWriter
        %F=getframe;writeVideo(writerObj,F);
    end
    %% plot trajectory
    if(opt.output.plotTrj)
        if(~ishandle(trjAxis)) % then the axis was probably deleted
            errordlg('SM_runsimulation: trajectory plot window not found (closed?). Terminating simulation.')
            break
        end
        axes(trjAxis)
        cla
        hold on
        
        for k=1:length(traj)
            if(~isempty(traj{k}))
                %plot3(trjAxis,traj{k}(:,3),traj{k}(:,4),traj{k}(:,5),'.b')
                ind=find((traj{k}(:,1)>=tFrame0).*(traj{k}(:,1)<=tFrame1));
                plot3(trjAxis,traj{k}(ind,3),traj{k}(ind,4),traj{k}(ind,5),'k-')
            end
        end
        clear k
        if(~isempty(xEmitted))
            plot3(trjAxis,xEmitted(:,1),xEmitted(:,2),xEmitted(:,3),'b.','markersize',1)
        end
        
        axis equal
        axis([[0 opt.camera.xrange_px+1 0 opt.camera.yrange_px+1]*opt.camera.pixLength opt.output.plotTrjZRange])
        %ZD=get(get(gca,'childr'),'zdata');ZD=[ZD{:}];ZD=[min(ZD) max(ZD)];
        %axis([[0 opt.camera.xrange_px+1 0 opt.camera.yrange_px+1]*opt.camera.pixLength ZD])
        %set(gca,'xlim',[0 opt.camera.xrange_px+1]*opt.camera.pixLength,'ylim',[0 opt.camera.yrange_px+1]*opt.camera.pixLength);
        
        view([20 10])         
        
        box on
        xlabel('x [nm]')
        ylabel('y [nm]')
        zlabel('z [nm]')
    end  
    %% pause and draw
    if(opt.output.plotTifMovie || opt.output.plotTrj)
        drawnow
        %pause(0.01)
%        if(nFrame==21)
%            makeSimImage
%            keyboard
%        end
    end
    %% save results and tracking truth
    if( mod(nFrame,opt.output.movieLength)==0 || nFrame==opt.output.maxFrames)
        tr=assembleTrajectories(tru,opt);
        tr.opt=opt;                 % also save opt structure and movie names
        tr.movieName=result_movieName;
        try
            save(fullResultName,'-struct','tr');
            disp(['saved trajectory truths to ' opt.output.resultFile ' at frame ' int2str(nFrame)])
        catch me
            display(['Something went wrong when saving to ' ...
                fullResultName ' saving in SMerror.mat instead'])
            save SMerror.mat
            throw(me)
        end
        clear savePath
        %%% debugging: save workspace as well
    end
end
%% save the truth one final time
tr=assembleTrajectories(tru,opt);
tr.opt=opt;                 % also save opt structure and movie names
tr.movieName=result_movieName;
try
    save(fullResultName,'-struct','tr');
    disp(['saved trajectory truths to ' opt.output.resultFile ' at frame ' int2str(nFrame)])
catch me
    display(['Something went wrong when saving to ' ...
        fullResultName ' saving in SMerror.mat instead'])
    save SMerror.mat
    throw(me)
end
clear savePath
%% close simulations indata files
fclose(fiT);
fclose(fiR);
end
%% convert collected ground truth into separate trajectories
function trTruth=assembleTrajectories(mvTruth,opt)
t0=toc;
fname=fieldnames(mvTruth); % field names to sort out

trTruth=struct;
for n=1:length(fname)
    trTruth.(fname{n})=cell(0);
end
tr_tStart=[];
    
    % transfer trajectories of length > 1
    % truTraj{i}      =[x(t) y(t) z(t) frameNr movieNr] %
    % em/detAverage{i}=[x(t) y(t) z(t) frameNr movieNr]
    % z(t)=0 for the detectionAverage
    % spotStats{i}=[nPh covXX covXY covYY];
    for k=1:length(mvTruth.positionSnapshots)
        if(~isempty(mvTruth.positionSnapshots{k}) && size(mvTruth.positionSnapshots{k},2)>1)
            tr_tStart(end+1)=mvTruth.detectionAverage{k}(1,4);
            for n=1:length(fname)
                trTruth.(fname{n}){end+1}   = mvTruth.(fname{n}){k};
            end
            
            % break up trajectories with multiple move numbers
            mv=trTruth.detectionAverage{end}(:,5);
            ind=find(diff(mv)~=0);
            if(~isempty(ind))
                kLong=length(trTruth.positionSnapshots);
                s12=[[1;ind+1] [ind;length(mv)]]; % each row is start and end of one movie period
                for r=2:size(s12,1) % create new trajectories
                    tr_tStart(end+1)=trTruth.detectionAverage{kLong}(s12(r,1),4);
                    for n=1:length(fname)
                        trTruth.(fname{n}){end+1}   = trTruth.(fname{n}){kLong}(s12(r,1):s12(r,2),:);
                    end
                end % only keep the first movie in trajectory kLong
                for n=1:length(fname)
                    trTruth.(fname{n}){kLong}   = trTruth.(fname{n}){kLong}(s12(1,1):s12(1,2),:);
                end
            end
        end
    end
    
    % sort on starting time
    [~,ind]=sort(tr_tStart);
    clear k tr_tStart
    
    % offset frame numbers
    % positionSnapshots{i}  =[x(t) y(t) z(t) frameNr movieNr] %
    % em/detAverage{i}      =[x(t) y(t) z(t) frameNr movieNr]
    % z(t)=0 for the detectionAverage
    % spotStats{i}=[nPh covXX covXY covYY];    
    for k=1:length(trTruth.positionSnapshots)
        mv=trTruth.detectionAverage{k}(:,5);
        fr=trTruth.detectionAverage{k}(:,4);
        if(opt.output.movieLength<Inf) % no point writing lots of Infs
            fr=fr-(mv-1)*opt.output.movieLength;
        end
        trTruth.positionSnapshots{k}(:,4)=fr;
        trTruth.emissionAverage{k}(:,4)  =fr;
        trTruth.detectionAverage{k}(:,4) =fr;        
    end
    clear k mv fr
    
    for n=1:length(fname)
        trTruth.(fname{n})=trTruth.(fname{n})(ind);
    end    
    trTruth.opt=opt; % save simulation options for future ref
    trTruth.runID=opt.runVariables.runID; % save generation time etc
    
    t1=toc;
    disp(['saving the truth: ' num2str(t1-t0) ' s.'])
end
