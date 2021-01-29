%%% Programa el metodo de cuadratura Gauss-Legendre para calcular la 
%%% funcion Bessel

[x,w] = gauleg(0,2*pi,100);
f =@(m,R) 1i^(-m)/(2*pi)*exp(1i*m*x).*exp(1i*R*cos(x));

%%% Calcula J0(8), J1(5) y J2(2) con 6 decimales de exactitud.
J = zeros(1,3);
J(1) = f(0,8)*w';
J(2) = f(1,5)*w';
J(3) = f(2,2)*w';

%%% Compara tu resultado con la rutina besselj
J2 = [besselj(0,8),besselj(1,5),besselj(2,2)];

format long
disp(J')
disp(J2')

%%% Escribe una rutina para graficar la funcion J1(2x) usando el metodo
%%% de cuadratura gauleg. Grafica la funcion en el intervalo de 
%%% x = [-10; 10]:
close all; set(0,'defaultTextInterpreter','latex');
n=1;
R = linspace(-10,10,500);
h = R(2)-R(1);
Y = zeros(1,length(R));
for r = -10:h:10
Y(n) = f(1,2*r)*w';
n = n +1;
end

figure(1); plot(R,real(Y),'.r')
ejes = gca;
ejes.FontSize = 10;
title('$J_1(2x)$ usando cuadratura gauleg','FontSize',15)
xlabel('$x$','FontSize',15);
ylabel('$J_1(2x)$','FontSize',15);
axis square;

