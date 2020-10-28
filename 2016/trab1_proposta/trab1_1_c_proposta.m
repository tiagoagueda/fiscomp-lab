%
% Fisica Computacional (2016-2017)
%
% Trabalho Pratico 1
% Problema 1.1 c) (Proposta de Resolucao)
%
% Author : Tiago A. Oliveira (tiago.agueda@ua.pt)
% Revisions :
% 2017/02/17 - File created.
%

clear all % clear all variable
close all % close all windows
clc       % clear terminal output

% universal constants:
g = 9.8;        % (m/s^2)   gravity acelaration

% fixed parameters:
v0 = 16.0;       % (m/s)     initial velocity
vlim = 6.8;     % (m/s)     terminal velocity
t0 = 0.0;       % (s)       initial time
tf = 2.0;       % (s)       final time
z0 = 1.0;       % (m)       starting position

% Euler Method - parameters:
h=0.1;          % $$h=t_{k+1}-t_{k}$$

% variable allocation
n = tf/h;                 % find size of matrix/vector
t = t0:h:tf;                % allocate and initialize t
vz = zeros(1,n);            % allocate v_z
vz(1) = v0;                 % initialize v_z
z = zeros(1,n);             % allocate z
z(1) = z0;                  % initialize z

% analytical solution
% vz_a=-vlim*tanh(g.*t/vlim);

for i=1:n
    vz(i+1)=vz(i)+h*g*(-1-abs(vz(i))*vz(i)/(vlim*vlim));
    z(i+1)=z(i)+vz(i)*h;
    if(z(i+1)<0)
        break
    end
end

% remove ununcessary elements
t=t(1:i+1);
z=z(1:i+1);
vz=vz(1:i+1);

t_solo=interp1(z(end-1:end),t(end-1:end),0,'linear');
vz_solo=interp1(z(end-1:end),vz(end-1:end),0,'linear');
fprintf('Collision with ground at t = %f s\n\n', t_solo)

% Correct final values
t(end)=t_solo;
z(end)=0;
vz(end)=vz_solo;


% ploting results
plot(t,z)
xlabel('time (s)')
ylabel('height (m)')

% Esta figura nao e pedida:
figure(2)
plot(t,vz)
xlabel('time (s)')
ylabel('velocity (m/s)')