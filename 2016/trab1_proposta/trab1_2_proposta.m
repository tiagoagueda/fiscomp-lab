%
% Fisica Computacional (2016-2017)
%
% Trabalho Pratico 1
% Problema 1.2 a)b)c) (Proposta de Resolucao)
%
% Author : Tiago A. Oliveira (tiago.agueda@ua.pt)
% Revisions :
% 2017/03/04 - Minor improvements.
% 2017/02/20 - File created.
%

clear all % clear all variable
close all % close all windows
clc       % clear terminal output

% universal constants:
g = [0 0 -9.8]; % (m/s^2)   gravity acelaration vector

% fixed parameters:
w_rpm = 0;      % a) No rotation!
%w_rpm = 3000;   % b) topspin 3000 RPM
%w_rpm = -3000;  % c) backspin 3000 RPM

z0 = 0.7;       % (m)       initial position (z)
v0 = 20.0;      % (m/s)     initial velocity
vangle = 5;     % (degrees) initial velocity angle with horizon
t0 = 0.0;       % (s)       initial time
tf = 10.0;      % (s)       final time
R = 67/2;       % (mm)      ball radius
m = 57;         % (g)       ball mass
rho = 1.225;    % (Kg/m^3)  air density

% unit conversion:
R = R*1e-3;     % mm -> m
m = m*1e-3;     % g -> Kg

% other parameters' calculation:
A = pi*R^2;             % (m^2)     section area
B = 0.5*A*rho;          %
wy=w_rpm*2*pi/60;       % (rad/s)   angular velocity (y component)
w=abs(wy);              % (rad/s)   angular velocity
w_vers=[0 sign(wy) 0];  % angular velocity vector
                        % The sign() function prevents zero division!


% Euler Method - parameters:
h=0.001;          % $$h=t_{k+1}-t_{k}$$

% variable allocation
t = t0:h:tf;                % allocate and initialize t
n = length(t);              % find size of vectors
r = zeros(n,3);             % position vector (x,y,z) - allocation
r(1,3) = z0;                % position vector initialization
v_vect = zeros(n,3);        % velocity vector (vx,vy,vz) - allocation 
v_vect(1,:) = [v0*cosd(vangle) 0 v0*sind(vangle)];    % - initialization
v = zeros(1,n);             % allocate v (modulus)
v(1) = v0;                  % initialize v (modulus)

for i=1:n-1
    S = R*w/v(i);
    CL = 1/(2.022+0.981/S);
    CD = 0.508 + (22.503+4.196*S^(-2.5))^(-0.4);
    
    a=g+(-B*CD*v(i)*v_vect(i,:)+B*CL*v(i)*cross(w_vers,v_vect(i,:)))/m;
    
    v_vect(i+1,:)=v_vect(i,:)+a*h;
    v(i+1)=norm(v_vect(i+1,:));
    r(i+1,:)=r(i,:)+v_vect(i,:)*h;
    if r(i+1,3)<0
        break
    end
end

% clean unwanted zeros
t = t(1:i+1);
r = r(1:i+1,:);
v_vect = v_vect(1:i+1,:);
v = v(1:i+1);

% Find maximum range
x_range=interp1(r(end-1:end,3),r(end-1:end,1),0,'linear');
t_range=interp1(r(end-1:end,3),t(end-1:end),0,'linear');
v_range=interp1(r(end-1:end,3),v(end-1:end),0,'linear');
fprintf('Rotation: %d rpm\n', w_rpm)
fprintf('Time of flight = %f s\n', t_range)
fprintf('Range = %f m\n', x_range)
fprintf('Impact velocity = %f m/s\n', v_range)

% Maximum heigh
% Find maximum value and its index
[z_max,ind]=max(r(:,3)); % We can use the value z_max.
% mas eu vou complicar e usar o valor maximo 
% de uma parabola que passa pelo valor maximo do vetor
% e pelos pontos anterior e posterior.
% O ficheiro maximo.m tem que estar no PATH.
aux=maximo(r(ind-1:ind+1,1),r(ind-1:ind+1,3));
z_max=aux(2);
fprintf('Maximum heigh = %f m\n\n', z_max)

% Refine same values for plotting
t(end)=t_range;
r(end,1)=x_range;
r(end,3)=0;
v(end)=v_range;

% ploting results
figure(1)
plot3(r(:,1),r(:,2),r(:,3))
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Trajectory')