%%%%%%%%% CUADRATURA GAUSSIANA %%%%%%%%%
a = 0; b = 6; n = 2^8;
[r1,w] = gauleg(a,b,n);
[R1,R2] = meshgrid(r1);

%%%%%%%%% PARAMETROS INICIALES %%%%%%%%%
lambda = 10.6*1e-3;
theta0 = 12.22e-3;
L = 409.2;                  %Distancia entre los elementos                  %Radio de curvatura
k = 2*pi/lambda;            %Número de onda
p = 1; m=1;

for n = 1:6

RR = [200*L 100*L 50*L 25*L 10*L 5*L];
R = RR(n);

x1 = (L-R)*theta0;
Q = L/R;
A = 1 - 2*Q - 2*(theta0./x1)*2*L*(1-Q);
B = 2*L*(1-Q);
%C =@(x1) (-2*theta0./x1)*(1-2*Q)-2/R; 
D = 1 - 2*Q;

%%%%%%%%% INTEGRAL %%%%%%%%%
H = @(l) k./B.*(-1i)^(l+1).*besselj(l,k.*R1.*R2./B).*...
    exp(1i*k./(2*B).*(A.*R1.^2+D.*R2.^2)).*R1.*exp(-1i*2*k*theta0*R1);
W = diag(w);
 K = H(0)*W;
 [V,Y] = eig(K);
 norm = max(abs(V(:,p)));
    
%%%%%%%%% PLOTS %%%%%%%%%
 
    
    figure(1);sgtitle('Eigenmodos para distintos valores de R para l=0')
    subplot(2,3,m); plot(r1,abs(V(:,p))/norm,'r')
    title('R = -' + string(RR(n)/L) +'L' )
    xlabel('r[mm]');
    ylabel('|U|');
    
    m = m +1;
end