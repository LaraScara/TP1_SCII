%Caso de estudio 2 - Motor CC apartados 1 y 2
clc, close all;

%Constantes
Laa = 366e-6;
J = 5e-9;
Ra = 55.6;
B = 0;
Ki = 6.49e-3;
Km = 6.53e-3;

X = [0; 0; 0; 0]; %valores iniciales
ii = 0;
t_etapa = 1e-7;
wRef = 2;
tF = 5; % 5 segundos
Tl = 0;


for t=0:t_etapa:tF
 if t==2        %retardo 2 seg para el torque
    Tl = 0.0014 % Con este torque se garantiza una velocidad nula
 end   
 ii = ii+1;
 X = modmotor(t_etapa, X, 12, Tl);
 x1(ii) = X(1); %w
 x2(ii) = X(2); %wp
 x3(ii) = X(3); %corriente
 acc(ii) = 12;
end

%Gráficas
t=0:t_etapa:tF;
figure
plot(t,x1);
title('Salida del sistema (velocidad angular)');
xlabel('Tiempo (segundos)');
ylabel('Velocidad angular (rad/s)')
legend('wr(t)','Location','southeast');
figure
plot(t,acc);
title('Entrada del sistema (tensión de armadura)');
xlabel('Tiempo (segundos)');
ylabel('Voltaje (volts)')
legend('va(t)','Location','southeast');
figure
plot(t,x3);
title('Corriente de armadura');
ylabel('Corriente (ampere)')
xlabel('Tiempo (segundos)');
legend('ia(t)','Location','southeast');
disp('Terminado')

% Método de Euler
function [X]=modmotor(t_etapa, xant, accion, Tl)
Laa=366e-6; J=5e-9; Ra=55.6; B=0; Ki=6.49e-3; Km=6.53e-3;
Va = accion;
h = t_etapa;
w = xant(1);
wp = xant(2);
ia = xant(3);
ip = xant(4);
ip = -Ra/Laa*ia - Km/Laa*w+Va/Laa;
ia = ia + h*ip;
wp = Ki/J * ia - B/J*w - Tl/J;
w = w + h*wp;
X = [w, wp, ia, ip];
end