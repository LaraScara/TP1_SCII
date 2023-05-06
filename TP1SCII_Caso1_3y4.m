%Caso de estudio 1 - Circuito RLC apartado 3 y 4
clc, clear all, close all;

%Importo gráficos de mediciones y grafico
mediciones=xlsread('Curvas_Medidas_RLC.xls');
%Columnas: Tiempo Corriente Vcap

%Gráfico de corriente
figure(1)
plot(mediciones(:,1),mediciones(:,2)); 
title('i_a variable de estado x_1')
xlabel('Tiempo [Seg.]');
ylabel('Corriente [Amp.]');

%Gráfico de voltaje de capacitor
figure(2)
plot(mediciones(:,1),mediciones(:,3)); 
title('V_c var. de est. x_2')
grid on;
xlabel('Tiempo [Seg.]');
ylabel('Voltaje [Volt]');

%Definicion del sistema RLC serie en variables de estado
%Variables de estado: x1=i, x2=vc

%Definicion de la entrada
paso=0.1/1000;
t=0:paso:(0.1-paso);
u=zeros(1,length(t));

signo=true;
for i=100:1:1000
    if mod(i,500)==0
       signo=not(signo);
    end
    if signo==1
        u(1,i)=12;
    end
    if signo==0
        u(1,i)=-12;
    end
end
figure
plot(t,u)

%Obtención de R L y C

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))

K=12;
t0=101;
t=50;

y1=mediciones(t0+t,3);
y2=mediciones(t0+2*t,3);
y3=mediciones(t0+3*t,3);

k1=y1/K-1; 
k2=y2/K-1; 
k3=y3/K-1; 
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1=-0.005/log(alfa1);
T2=-0.005/log(alfa2);
T3=beta*(T1-T2)+T1; %No hay cero, no se usa

s=tf('s');
G=12/((T1*s+1)*(T2*s+1));
[num,den]=tfdata(G,'v');

%Se sabe que Vc/vin (s) = (1/(L*C)) / (s^2+R*s/L+1/(L*C))
den_norm=den/(den(1))
num_norm=num/(12*den(1))
L=0.1;
R=L*den_norm(2);
Cap=1/(L*den_norm(3));

figure
[yaprox,taprox]=lsim(G,u/12,t);
plot(0.015,y1,'x');
hold on;
plot(0.02,y2,'x');
plot(0.025,y3,'x');
plot(taprox,yaprox);
title('Respuesta del voltaje del capacitor aproximada');
xlabel('Tiempo');
legend({'vc(t1)','vc(2t1)','vc(3t1)','vc(t)'},'Location','southeast')

%verifico que el sistema está bien aproximado
figure
plot(mediciones(:,1),mediciones(:,3));
hold on;
plot(taprox,yaprox); %Son muy parecidos, la aproximación es válida
legend({'vc(t) real','vc(t) aproximada'},'Location','southeast')
title('Contraste de respuestas');
xlabel('Tiempo');

%se verifica la corriente?

%Matrices
A=[-R/L -1/L; 1/Cap 0];
B=[1/L; 0];
C=[1; 0]';
D=0;

%Definicion de la ecuación de estado y de salida (salida de corriente)
sys1=ss(A,B,C,D)
figure
[yout,yt]=lsim(sys1,u,t);
plot(yt,yout);
hold on;
plot(mediciones(:,1),mediciones(:,2)); %se verifica con los componentes obtenidos
legend({'i(t) aproximada','i(t) real'},'Location','southeast')
title('Contraste de respuestas');
xlabel('Tiempo');