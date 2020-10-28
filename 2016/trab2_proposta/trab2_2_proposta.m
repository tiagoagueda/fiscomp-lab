clc
clear all
close all

% Valores iniciais 
x0=0.47;
y0=0;
vx0=0;
vy0=8.2;

% Parametros do metodo de Euler-Cromer
tfin=0.5; % 2 anos da Terra da' mais que oito periodos
h=0.0001;

% Constantes Fisicas
GMs=4*pi^2;

% Iniciar variaveis
t=0:h:tfin;
N=numel(t);
x=zeros(1,N);
x(1)=x0;
y=zeros(1,N);
y(1)=y0;
r=zeros(1,N);
r(1)=norm([x0 y0]);
vx=zeros(1,N);
vx(1)=vx0;
vy=zeros(1,N);
vy(1)=vy0;
v=zeros(1,N); % Nao e' necessario calcular v
v(1)=norm([vx0 vy0]);
ang=zeros(1,N);
ang(1)=mod(atan2(y0,x0),2*pi); % = 0

% Euler-Cromer

    for k=1:N-1
        vx(k+1)=vx(k)-GMs*x(k)/r(k)^3*h;
        vy(k+1)=vy(k)-GMs*y(k)/r(k)^3*h;
        v(k+1)=norm([vx(k+1) vy(k+1)]);
        x(k+1)=x(k)+vx(k+1)*h;
        y(k+1)=y(k)+vy(k+1)*h;
        r(k+1)=norm([x(k+1) y(k+1)]);
        ang(k+1)=mod(atan2(y(k+1),x(k+1)),2*pi);
        if ang(k+1)<ang(k)
            break
        end    
    end
    
    % Elimina os zeros finais
    N=k+1; % Novo numero de pontos
    t=t(1:N);
    x=x(1:N);
    y=y(1:N);
    r=r(1:N);
    vx=vx(1:N);
    vy=vy(1:N);
    v=v(1:N);
    ang=ang(1:N);
    ang(N)=ang(N)+2*pi; % Para as contas darem certo na alinea b), 
                        % o angulo nao pode dar um salto de -2 pi 
    
    % a)                    
    figure(1)
    plot(x,y)
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    axis([-0.5 0.5 -0.5 0.5])
    xlabel('x'),ylabel('y');
    title('Alinea a)')
    
    % b)
    T=interp1(ang(end-1:end),t(end-1:end),2*pi);
    fprintf('Periodo = %d\n\n', T)
    
    % c)
    Nmeio=round(N/2);
    Delta_Area=zeros(1,Nmeio);
    for k=1:Nmeio
        Delta_Area(k)=((r(k)+r(k+1))/2)^2*(ang(k+1)-ang(k));
    end    
    t=t(1:Nmeio);
    figure(2)
    plot(t,Delta_Area)
    xlabel('t'),ylabel('Variacao de Area');
        title('Alinea c)')

    
