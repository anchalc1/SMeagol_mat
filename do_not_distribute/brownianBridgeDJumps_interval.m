function xi = brownianBridgeDJumps_interval(x0,x1,a,b,ti,dt,D,tD)
% brownian bridge interpolation from x0 to x1, at time 0 <= ti <= dt, with
% a pice-wise constant diffusion constant.  xi is the interpolated
% position, which is restricted to lie in the interval a < xi < b.
%
% D(t) is specified by a list of diffusion constants and jump times tD,
% with the convention that D=D(m) when tD(m) <= t < tD(m+1) (hence
% tD(1)=0 is required). Omitting tD is equivalent to specifying tD=0, i.e., 
% the situation with a single diffusion constant.

% ML 2013-10-22

% parameter check
if(~exist('tD','var') || isempty(tD)), tD=0; end
if(length(tD) ~= length(D))
    error('brownianBridgeDJump_interval: D and tD must have equal lengths'); 
end
if(tD(1) ~= 0)
    error('brownianBridgeDJump_interval: tD(1)=0 required'); 
end
if( a ~< b)
	error('brownianBridgeDJump_interval: a < b required'); 
end

% construct cumulative variance Vt
Vt_t=[tD dt]; 
Vt_var=[0 cumsum(2*D.*diff(Vt_t))]; 
Vt=interp1(Vt_t,Vt_var,ti,'linear');
V1=interp1(Vt_t,Vt_var,dt,'linear');

% simulate a diffusion for the requested interpolation time
mu_t =x0+(x1-x0)*Vt/V1;
sig_t=sqrt(Vt*(1-Vt/V1));
u=rand;

xi = mu_t + sig_t*erfinv( (1-u)*erf((a-mu_t)/sig_t)+u*erf((b-mu_t)/sig_t));

end


