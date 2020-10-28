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

n = input(['Escolha o metodo:\n  1-Euler,\n  2-Euler-Cromer,\n  '...
    '3-Euler implicito sem linsolve,\n  4-Euler implicito com linsolve,\n  '...
    '5-Crank-Nicolson com linsolve.\n']);

switch n
    
    case 1
        fprintf('\nEuler\n\n')
        for k=1:N-1
            a=-K*x(k)/m;
            vx(k+1)=vx(k)+a*h;
            x(k+1)=x(k)+vx(k)*h;
        end
        
    case 2
        fprintf('\nEuler-Cromer\n\n')
        for k=1:N-1
            a=-K*x(k)/m;
            vx(k+1)=vx(k)+a*h;
            x(k+1)=x(k)+vx(k+1)*h;
        end
        
    case 3
        fprintf('\nEuler implicito sem linsolve\n\n')
        aux=1+w2*h^2;
        for k=1:N-1
            x(k+1)=(x(k)+vx(k)*h)/aux;
            vx(k+1)=vx(k)-w2*x(k+1)*h;
        end
        
    case 4
        fprintf('\nEuler implicito com linsolve\n\n')
        A=[1 -h; w^2*h 1];
        for k=1:N-1
            b=[x(k); vx(k)];
            aux=linsolve(A,b); 
            x(k+1)=aux(1); 
            vx(k+1)=aux(2); 
        end
        
     case 5
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

E=0.5*(K*x.^2+m*vx.^2);
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


