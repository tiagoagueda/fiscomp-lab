%
% Fisica Computacional (2016-2017)
%
% Trabalho Pratico 1
% Problema 1.1 a) (Proposta de Resolucao)
%
% Author : Tiago A. Oliveira (tiago.agueda@ua.pt)
% Revisions :
% 2017/02/16 - Updated revision.
% 2017/02/04 - File created.
%

clear all % clear all variable
close all % close all windows
clc       % clear terminal output

% universal constants:
g = 9.8;        % (m/s^2)   gravity acelaration

% fixed parameters:
v0 = 0.0;       % (m/s)     initial velocity
vlim = 6.8;     % (m/s)     terminal velocity
t0 = 0.0;       % (s)       initial time
tf = 2.0;       % (s)       final time

% Euler Method - parameters:
h=0.1;          % $$h=t_{k+1}-t_{k}$$

% variable allocation
n = tf/h;                 % find size of matrix/vector
t = t0:h:tf;                % allocate and initialize t
vz = zeros(1,n);            % allocate v_z
vz(1) = v0;                 % initialize v_z

% analytical solution
vz_a=-vlim*tanh(g.*t/vlim);

for i=1:n
    vz(i+1)=vz(i)+h*g*(-1-abs(vz(i))*vz(i)/(vlim*vlim));
end

% ploting results
plot(t,vz,'r',t,vz_a,'b')
legend('numerical','analytical','Location','Best')
xlabel('time (s)')
ylabel('velocity (m/s)')