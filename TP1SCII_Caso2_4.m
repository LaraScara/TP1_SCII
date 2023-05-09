%Caso de estudio 2 - Motor CC apartado 4
clc, close all;

%Valores iniciales y referencias
X=[0; 0; 0; 0; 0]; 
titaRef=1;
Tl=0;
t_etapa=1e-7;
tF=3;

%Constantes del PID en tiempo continuo
Kp=50;
Ki=300;
Kd=3;

%Conversion a PID en tiempo discerto
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(tF/t_etapa,1);

u=0;
ii=0;
for t=0:t_etapa:tF
 if t==0.5
    Tl=0.0007;
 end   
 ii=ii+1;
 k=ii+2;
 X=modmotor(t_etapa, X, u, Tl);
 e(k)=titaRef-X(5); %ERROR
 u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
 if u>12
     u=12;
 end
 if u<-12
     u=-12;
 end    
 x1(ii)=X(1);%w
 x2(ii)=X(2);%wp
 x3(ii)=X(3);%corriente
 x4(ii)=X(4);%corriente punto
 x5(k)=X(5);%tita, posición
 acc(k)=u;
end


t=0:t_etapa:tF+2*t_etapa;
figure
plot(t,e);
title('Error');
xlabel('Tiempo (segundos)');
legend('e(t)','Location','southeast');
axis([0 2.5 -0.2 1.2]);
figure
plot(t,acc);
title('Entrada del sistema (tensión de armadura)');
xlabel('Tiempo (segundos)');
ylabel('Voltaje (volts)');
legend('va(t)','Location','southeast');
axis([0 2.5 -15 15]);
figure
plot(t,x5);
title('Posición angular del motor');
xlabel('Tiempo (segundos)');
ylabel('Posición angular (rad)');
legend('θ(t)','Location','southeast');
axis([0 2.5 -0.2 1.2]);

disp('Terminado');

function [X]=modmotor(t_etapa, xant, accion, Tl)
Laa=366e-6; J=5e-9;Ra=55.6;B=0;Ki=6.49e-3;Km=6.53e-3;
Va=accion;
h=t_etapa;
w= xant(1);
wp = xant(2);
ia = xant(3);
ip = xant(4);
tita = xant(5);
ip = -Ra/Laa*ia - Km/Laa*w+Va/Laa;
ia = ia + h*ip;
wp = Ki/J * ia - B/J*w - Tl/J;
w = w + h*wp;
tita = tita + h*w;
X=[w, wp, ia, ip, tita];
end
