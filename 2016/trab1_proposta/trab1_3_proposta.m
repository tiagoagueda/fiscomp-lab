clc
clear all
close all

% Lancamento segundo xx, zz vertical

% Input que pode ser modificado
wz_rpm=+600;
v0=80/3.6;
ang0=10;
% Calculos auxiliares
wz = wz_rpm*2*pi/60;
w_vec = [0 0 wz];

% Metodo de Euler
h=0.001;
tfin=2;
% tfin nao tem que ser multiplo de h (t nao chega a tfin)

% Constantes do problema
g=[0 0 -9.8];
rho=1.225; 
m=0.450; 
R=0.70/2/pi; 
A=pi*R^2; 
CM=1; 
B=CM*rho*A*R/2;

% Iniciar Variaveis
t=0:h:tfin;
N=numel(t);
r_vec=zeros(N,3); % Vetor posicao
v_vec=zeros(N,3);
v_vec(1,:)=[v0*cosd(ang0) 0 v0*sind(ang0)];


for k=1:N-1
    v=norm(v_vec(k,:));
    if v<= 9
        aux=0.015*v^2;
    elseif v>20
        aux=-4.025+0.323*v;
    else
        aux=0.25147+0.17431*v-0.01384*v^2+0.00054*v^3;
    end
    
    FD=-aux/v*v_vec(k,:);
    FL=B*cross(w_vec,v_vec(k,:));
    
    a=(FD+FL)/m+g;

    v_vec(k+1,:)=v_vec(k,:)+a*h;
    r_vec(k+1,:)=r_vec(k,:)+v_vec(k,:)*h;
    
    if r_vec(k+1,3)<0 
        break
    end    
end

% Elimina os zeros finais
t=t(1:k+1);
r_vec=r_vec(1:k+1,:);
v_vec=v_vec(1:k+1,:);

% Alcance
t_alc=interp1(r_vec(k:k+1,3),t(k:k+1),0,'linear');
v_alc=interp1(r_vec(k:k+1,3),v_vec(k:k+1,:),0,'linear');
r_alc=interp1(r_vec(k:k+1,3),r_vec(k:k+1,:),0,'linear');
fprintf('Distancia percorrida segundo x = %f m\n', r_alc(1))
fprintf('Desvio lateral para a esquerda = %f m\n\n', r_alc(2))

% Acerta alguns pontos finais para graficos
t(end)=t_alc;
r_vec(end,:)=r_alc;
v_vec(end,:)=v_alc;

figure(1)
subplot(2,1,2)
plot(r_vec(:,1),r_vec(:,2))
xlabel('x');ylabel('y')
subplot(2,2,1)
plot(r_vec(:,1),r_vec(:,3))
xlabel('x');ylabel('z')
subplot(2,2,2)
plot(t,r_vec(:,3))
xlabel('t');ylabel('z')

figure(2)
plot3(r_vec(:,1),r_vec(:,2),r_vec(:,3))
xlabel('x');ylabel('y');zlabel('z')