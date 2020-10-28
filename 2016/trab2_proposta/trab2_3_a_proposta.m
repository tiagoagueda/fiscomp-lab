clc
clear all
close all


% Input que pode ser modificado
alfa=-0.1;
x0=1.0;
vx0=1.0;

% Parametros dos metodos numericos
h=0.01;
tfin=50;

% Constantes do problema
K=1; % K maiusculo!
m=1;
% Calculos auxiliares
w=sqrt(K/m);
w2=K/m;

% Iniciar Variaveis
t=0:h:tfin;
N=numel(t);
x=zeros(N,1);
x(1)=x0;
vx=zeros(N,1);
vx(1)=vx0;

% Euler-Cromer
for k=1:N-1
    a=-K*(x(k)+2*alfa*x(k)^3)/m;
    vx(k+1)=vx(k)+a*h;
    x(k+1)=x(k)+vx(k+1)*h;
end
E=0.5*(K*(x.^2+alfa*x.^4)+m*vx.^2); % Nao e' pedido


%Localiza maximos
imax=0;
for k=2:N-1
    if and(x(k+1)-x(k)<=0,x(k)-x(k-1)>=0)
        imax=imax+1;
        aux=maximo(t(k-1:k+1),x(k-1:k+1));
        tmax(imax)=aux(1);
        xmax(imax)=aux(2);
    end
end
nmax=imax; % Desnecessario
plf=polyfit(1:nmax,tmax,1);
T=plf(1);
A=mean(xmax);

% Figuras
figure(1)
subplot(2,2,1)
plot(t,x)
xlabel('\it t');ylabel('\it x')
subplot(2,2,2)
plot(t,vx)
xlabel('\it t');ylabel('\it vx')
subplot(2,2,3)
plot(x,vx)
xlabel('\it x');ylabel('\it vx')
subplot(2,2,4)
plot(t,E)
xlabel('\it t');ylabel('\it E')


% Mostra resultados
fprintf('Período metodo. numerico: %d \n\n', T)
fprintf('Amplitude metodo. numerico: %d \n\n', A)
