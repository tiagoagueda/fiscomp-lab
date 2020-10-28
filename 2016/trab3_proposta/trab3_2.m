%
% Fisica Computacional (2016-2017)
%
% Trabalho Pratico 3
% Problema 3.1 (Proposta de Resolucao)
%
% Author : Tiago A. Oliveira (tiago.agueda@ua.pt)
% Revisions :
% 2017/03/09 - File created.
%

clear all % clear all variable
close all % close all windows
clc       % clear terminal output

x0 = 1.0;   % (m)   - initial position
vx0 = 0.0;  % (m/s) - initial velocity
K = 16;     % (N/m) - elastic constant
m = 1.0;    % (Kg)  - pendulum mass

% methods related parameters
t0 = 0.0;   % (s)   - simulation initial time
tf = 10.0;  % (s)   - simulation final time
h = 0.10;   % (s)   - time increment

% variables allocation and initialization
t=t0:h:tf;      % allocate and initializate time vector
N=length(t);    % find size of vectors
% - for Runge-Kutta solutions
xRG=zeros(N,1);   % allocate position (x) vector
xRG(1)=x0;        % initialize position (x) vector
vxRG=zeros(N,1);  % allocate velocity (vx) vector
vxRG(1)=vx0;      % initialize velocity (vx) vector
% - for Euler solutions
xE=xRG;
vxE=vxRG;

% Functions definition:
fx = @(V) V;      % dx/dt=v
fv = @(X) -K*X/m; % dv/dt=a=-K*x/m


% Methods cycle
for i=1:N-1
    %Runge-Kutta
    r1v=fv(xRG(i));         % Check 'Aula Teorica 3'
    r1x=fx(vxRG(i));         % Slide 8
    r2v=fv(xRG(i)+r1x*h/2); % -
    r2x=fx(vxRG(i)+r1v*h/2); % -
    vxRG(i+1)=vxRG(i)+r2v*h;  % Check 'Aula Teorica 3'
    xRG(i+1)=xRG(i)+r2x*h;  % Slide 9

    % Euler
    vxE(i+1)=vxE(i)-K*xE(i)*h/m;
    xE(i+1)=xE(i)+vxE(i)*h;
end

% Mechanical Energy
ERG=0.5*K*xRG.^2+0.5*m*vxRG.^2;  % Runge-Kutta
EE=0.5*K*xE.^2+0.5*m*vxE.^2;     % Euler

% Analitical Solution
w=sqrt(K/m);
Eas=0.5*(K*x0^2+m*vx0^2)*ones(N,1);
xas=x0*cos(w*t);
vas=-w*x0*sin(w*t);


plot(t,xRG,'-k',t,xE,'-r')