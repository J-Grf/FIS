%
% - Amp factor

close all; clear all;

Nv=[ 10 50  200 500  1000 ];
% Nv=[400];
mu=zeros(1,length(Nv));
eps = 0.0001;

for l=1:length(Nv)
   
N=Nv(l);
dx   = 2.0*pi/(N-1);
dy   = dx;

theta1    = -pi:dx:pi;
theta2    = -pi:dy:pi;

[T1,T2]=meshgrid(theta1,theta2);


% -- "isotropic case"
am  = exp(1i*T1) +exp(1i*T2);
ctr = 4.0 - exp(-1i*T1) -exp(-1i*T2);

% -- "anisotropic case"
% am  = eps*exp(1i*T1)+exp(1i*T2);
% ctr = (2.0 +2.0*eps)- eps*exp(-1i*T1) -exp(-1i*T2);

% -- "anisotropic case with line relaxation"
% am  = eps*exp(1i*T1);
% ctr = (2.0 +2.0*eps)- eps*exp(-1i*T1) -exp(-1i*T2) -exp(1i*T2);

S = am./ctr;

surf(T1,T2,abs(S))
set(gca,'FontSize',18); pause(1)

mu(l)=max(max(abs(S(abs(T1)>pi/2.0 | abs(T2)>pi/2.0))));

end

if length(Nv) > 1
    plot(Nv,mu,'--ks')
    set(gca,'FontSize',18);
    xlabel('N'); ylabel('max |\lambda|');
    title('Maximum in high-frequency range')
end
