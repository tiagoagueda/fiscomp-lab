clc
clear all
close all


% Input que pode ser modificado
alfa=-0.1;
x0_vec=0.1:0.1:2.0;
NA=numel(x0_vec);
vx0=0.0; % Tem que ser 0 para que A=x0

% Parametros dos metodos numericos
h=0.01;
tfin=50;

% Constante0s do problema
K=1; % K maiusculo!
m=1;

% Iniciar Variaveis
t=0:h:tfin;
N=numel(t);

% Euler-Cromer para varias amplitudes
for iA=1:NA
    
    % Iniciar Variaveis
    x=zeros(N,1); % Nao era preciso, mas e' bom habito limpar
    x(1)=x0_vec(iA);
    vx=zeros(N,1); % Nao era preciso, mas e' bom habito limpar
    vx(1)=vx0;
    
    % Euler-Cromer
    for k=1:N-1
        a=-K*(x(k)+2*alfa*x(k)^3)/m;
        vx(k+1)=vx(k)+a*h;
        x(k+1)=x(k)+vx(k+1)*h;
    end
    
    %Localiza maximos
    imax=0;
    clear tmax xmax
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
    T(iA)=plf(1);
    A(iA)=mean(xmax);
end

% Figuras
figure(1)
plot(A,T,'*')
xlabel('\it A')
ylabel('\it T')




