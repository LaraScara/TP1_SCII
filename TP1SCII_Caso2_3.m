%Caso de estudio 2 - Motor CC apartado 3
clc, close all;

%Importo gráficos de mediciones y grafico
mediciones=xlsread('Curvas_Medidas_Motor_2023.xls');
%Columnas: Tiempo , velocidad angular , corriente , tension , torque

%Graficos de curvas con los datos extraidos de la tabla
figure
plot(mediciones(:,1),mediciones(:,2));%grafico de velocidad angular
title('Salida del sistema (velocidad angular)');
xlabel('Tiempo (segundos)');
ylabel('Velocidad angular (rad/s)');
legend('wr(t)','Location','northeast');
axis([0 0.6 -45 210]);

figure
plot(mediciones(:,1),mediciones(:,3)); %grafico de corriente de armadura
title('Corriente de armadura');
ylabel('Corriente (ampere)');
xlabel('Tiempo (segundos)');
legend('ia(t)','Location','northeast');
axis([0 0.6 -0.25 0.15]);

figure
plot(mediciones(:,1),mediciones(:,4)); %grafico de tension de entrada
title('Entrada del sistema (tensión de armadura)');
xlabel('Tiempo (segundos)');
ylabel('Voltaje (volts)');
legend('va(t)','Location','northeast');
axis([0 0.6 -15 15]);

figure
plot(mediciones(:,1),mediciones(:,5));%grafico del torque
title('Torque de carga');
xlabel('Tiempo (segundos)');
ylabel('Torque (Nm)');
legend('Tl(t)','Location','northeast');
axis([0 0.6 -0.0012 0.0001]);

%Obtención de FdT sin torque

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G1(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))
K = 198.248802246762;
a0=350;
a1=7000;

t0=mediciones(100,1);
t1=mediciones(a0,1);
t2=mediciones(a0+a1,1);
t3=mediciones(a0+2*a1,1);

y1=mediciones(a0,2);
y2=mediciones(a0+a1,2);
y3=mediciones(a0+2*a1,2);

k1=y1/K-1; 
k2=y2/K-1; 
k3=y3/K-1; 
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1=-(t1-t0)/log(alfa1)
T2=-(t1-t0)/log(alfa2)
T3=beta*(T1-T2)+T1;     %No hay cero, no se usa

s=tf('s');
G1=exp(-t0*s)*K/((T1*s+1)*(T2*s+1))

%Gráfica de los puntos seleccionados para Chen
figure
plot(mediciones(:,1),mediciones(:,2));
hold on;
plot(t1,y1,'o');
plot(t2,y2,'o');
plot(t3,y3,'o');
title('Salida del sistema (velocidad angular)');
grid on;
xlabel('Tiempo (segundos)');
ylabel('Velocidad angular (rad/s)');
legend('wr(t)','Location','northeast');
axis([0 0.6 -45 210]);
grid off;

%Obtención de FdT con torque

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G2(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))
K2 = -27;
b0=15500;
b1=1500;

tt=mediciones(15306,1);
t4=mediciones(b0,1);
t5=mediciones(b0+b1,1);
t6=mediciones(b0+2*b1,1);

y4=mediciones(b0,2);
y5=mediciones(b0+b1,2);
y6=mediciones(b0+2*b1,2);

k4=y4/K2 - 1; 
k5=y5/K2- 1;
k6=y6/K2 - 1;
b2=4*k4^3*k6-3*k4^2*k5^2-4*k5^3+k6^2+6*k4*k5*k6;
alfa3=(k4*k5+k6-sqrt(b2))/(2*(k4^2+k5));
alfa4=(k4*k5+k6+sqrt(b2))/(2*(k4^2+k5));
beta2=(2*k4^3+3*k4*k5+k6-sqrt(b2))/(sqrt(b2));
T4=-t4/log(alfa3)
T5=-t4/log(alfa4)
T6=beta2*(T4-T5)+T4

G2=exp(-tt*s)*(K-K2)+(T3*s+1)/((T4*s+1)*(T5*s+1))

%Gráfica de los puntos seleccionados para Chen
figure
plot(mediciones(:,1),mediciones(:,2));
hold on;
plot(t4,y4,'o');
plot(t5,y5,'o');
plot(t6,y6,'o');
title('Salida del sistema (velocidad angular)');
grid on;
xlabel('Tiempo (segundos)');
ylabel('Velocidad angular (rad/s)');
legend('wr(t)','Location','northeast');
axis([0 0.6 -45 210]);
grid off;

G=G1-G2

%Comparación de los valores aproximados y los valores reales
figure
plot(mediciones(:,1),mediciones(:,2));
hold on
step(G,'r'), title('Contraste de respuestas');
xlabel('Tiempo (segundos)');
ylabel('Velocidad angular (rad/s)');
legend({'wr(t) real','wr(t) aproximada'},'Location','northeast');
axis([0 0.6 -45 210]);

disp('Terminado')
