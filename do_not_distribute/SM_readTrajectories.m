% traj=SM_readTrajectories(RIF,fig)
%
% extract all single particle trajectories from the trajectory and reaction
% files associated with runinput file RIF. Output format is a cell vector
% for each particle, with each row containing [t id s x y z]. 
%
% If the optional argument fig is given, the trajectories are plotted
% in a 3D line plot n figure fig (or a new figure if fig=[]).
%
% ML 2013-11-04

function traj=SM_readTrajectories(RIF,fig)
opt=SM_getOptions(RIF);

%% read trajectory information
fiT=fopen(opt.trj.positionFile,'r');
fiR=fopen(opt.trj.reactionFile,'r');

% from trajectory file
traj=cell(1);
while(true)    
    [t,pid,state,x]=SM_readLine(fiT,opt.trj.A,opt.trj.b);
    if(isempty(t))
        break
    end
    for k=1:length(pid)
        if(pid(k)>length(traj) || isempty(traj{k}))            
            traj{pid(k)}=[t state(k) x(k,:)];
        else
            traj{pid(k)}(end+1,:)=[t state(k) x(k,:)];
        end
    end
end

% read in reactions
while(true)
    [t,pid,state,x]=SM_readLine(fiR,opt.trj.A,opt.trj.b);
    if(isempty(t))
        break
    end
    traj{pid}(end+1,:)=[t state x];
end

for k=1:length(traj)
   [~,ind]=sort(traj{k}(:,1));
   traj{k}=traj{k}(ind,:);    
end

fclose(fiT);
fclose(fiR);

%% plot?
if(exist('fig','var'))
    if(isempty(fig))
        figure
    else
        figure(fig)
    end
    
    clf
    hold on
    
    for k=1:length(traj)
        plot3(traj{k}(:,3),traj{k}(:,4),traj{k}(:,5))
    end
    
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
