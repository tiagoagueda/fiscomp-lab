%
% Fisica Computacional (2016-2017)
%
% Trabalho Pratico 1
% Problema 1.1 b) (Proposta de Resolucao)
%
% Author : Tiago A. Oliveira (tiago.agueda@ua.pt)
% Revisions :
% 2017/02/27 - File created.
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
tf = 0.5;       % (s)       final time

% Euler Method - parameters:
h=[0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002 0.0001];

% analytical solution
vz_a=-vlim*tanh(g.*tf/vlim);
    
for ih=1:length(h)
    % variable allocation
    n = tf/h(ih);                 % find size of matrix/vector
    t = t0:h(ih):tf;                % allocate and initialize t
    vz = zeros(1,n);            % allocate v_z
    vz(1) = v0;                 % initialize v_z



    for i=1:n
    vz(i+1)=vz(i)+h(ih)*g*(-1-abs(vz(i))*vz(i)/(vlim*vlim));
    end
    
    vz_fin = vz(end);
    erro(ih) = vz_fin-vz_a;
    log_mod_erro(ih)=log(abs(erro(ih)));
end

% ploting results
plot(log(h),log_mod_erro,'*')
%legend('numerical','analytical','Location','Best')
xlabel('log(h)')
ylabel('log(|error|)')

aux=polyfit(log(h),log_mod_erro,1);
fprintf('Exponent = %d \n\n', aux(1))