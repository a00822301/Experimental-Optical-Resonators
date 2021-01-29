%***********************************************************************
%
% El siguiente codigo es una traduccion de fortran a Matlab de la subrutina
% GAULEG, la cual se encuentra en el libro:
%     NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING
%
% Dados los limites x1,x2, y n la funcion regresa las abcisas y pesos del
% metodo de cuadratura de Gauss-Legendre para n puntos.
%
%                                      Raul Hernandez Aranda, 04/2020
%***********************************************************************
function [x,w] = gauleg(x1,x2,n)

EPS = 3e-14;

m = (n + 1)/2;      %Las raices son simetricas en el intervalo. Solo encontramos la mitad de ellas.
xm = 0.5*(x2 + x1);
xl = 0.5*(x2 - x1);

aux = 1;

for ii=1:m                              %Ciclo para las raices deseadas
    z = cos(pi*(ii - .25)/(n + .5));    %Se propone una raiz inicial
    z1 = z + 1;
    while abs(z - z1) > EPS             %Se hace el ciclo para obtener
        p1 = 1;                         %el polinomio de Legendre evaluado
        p2 = 0;                         %en z
        for jj=1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj - 1)*z*p2-(jj - 1)*p3)/jj;
        end
        %pp es el polinomio de Legendre deseado. Ahora se calcula su derivada
        %por medio de una relacion estandar que involucra a p2, el polinomio
        %de un orden menor.
        pp = n*(z*p1-p2)/(z*z-1);
        z1 = z;
        z = z1-p1/pp;                   %Metodo de Newton
    end
    x(ii) = xm - xl*z;                  %Se escala la raiz al intervalo deseado
    x(n+1-ii) = xm + xl*z;              %Se colocan las raices simetricas
    w(ii) = 2*xl/((1-z*z)*pp*pp);       %Se calculan los pesos
    w(n+1-ii) = w(ii);                  %Se colocan los pesos symetricos
end