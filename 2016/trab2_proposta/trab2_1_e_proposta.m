clc
clear all
close all

% Input que pode ser modificado
x0=1.0;
vx0=0.0;

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




n = input(['Escolha o metodo:\n  1-Euler-Cromer,\n  '...
        '2-Crank-Nicolson com linsolve.\n']);

switch n
    
    case 1

        fprintf('\nEuler-Cromer\n\n')
        for k=1:N-1
            a=-K*x(k)/m;
            vx(k+1)=vx(k)+a*h;
            x(k+1)=x(k)+vx(k+1)*h;
        end
        
     case 2
        fprintf('\nCrank-Nicolson com linsolve\n\n')
        A=[1 -h/2; w^2*h/2 1];
        for k=1:N-1
            b=[x(k)+vx(k)*h/2; vx(k)-w^2*x(k)*h/2];
            aux=linsolve(A,b); 
            x(k+1)=aux(1); 
            vx(k+1)=aux(2); 
        end

    otherwise
        fprintf('Nope')
        return
end


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

% Solução analitica
Asa=sqrt(x0^2+m*vx0^2/K);
Tsa=2*pi/w;
 
% Mostra resultados
fprintf('Período metod. numerico: %d \n', T)
fprintf('Período solu. analítica: %d \n', Tsa)
fprintf('Um a dividir pelo outro: %d \n\n', T/Tsa)
fprintf('Amplitude metd. numerico: %d \n', A)
fprintf('Amplitude sol. analítica: %d \n', Asa)
fprintf('Uma a dividir pelo outra: %d \n\n', A/Asa)