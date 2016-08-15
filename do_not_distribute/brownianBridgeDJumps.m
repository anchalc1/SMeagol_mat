% xi = brownianBridgeDJumps(x0,x1,ti,dt,D,tD)
%
% brownian bridge interpolation from x0 to x1, at times 0 <= ti <= dt, with
% a pice-wise constant diffusion constant.  xi are the interpolated
% positions.
%
% D(t) is specified by a list of diffusion constants and jump times tD,
% with the convention that D=D(m) when tD(m) <= t < tD(m+1) (hence
% tD(1)=0 is required). Omitting tD is equivalent to specifying tD=0, i.e., 
% the situation with a single diffusion constant.

% ML 2013-10-17

function xi = brownianBridgeDJumps(x0,x1,ti,dt,D,tD)

% parameter check
if(~exist('tD','var') || isempty(tD)), tD=0; end
if(length(tD) ~= length(D))
    error('brownianBridgeDJump: D and tD must have equal lengths'); 
end
if(tD(1) ~= 0)
    error('brownianBridgeDJump: tD(1)=0 required'); 
end

% construct cumulative variance Vt
Vt_t=[tD dt]; 
Vt_var=[0 cumsum(2*D.*diff(Vt_t))]; 
Vt_f=@(t)(interp1(Vt_t,Vt_var,t,'linear'));

% simulate a diffusion for the requested interpolation times
dx_std=sqrt(diff(Vt_f([0 ti(1:end) dt]))); % std of step lengths
dx=dx_std.*randn(size(dx_std));
x=cumsum(dx);

% impose end-point constraints
X0=x(1);
X1=x(end);
xi=x(1:end-1)+(x0-X0)+(x1-X1-x0+X0)*Vt_f(ti)/Vt_f(dt);
end

