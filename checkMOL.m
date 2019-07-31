% * * * * u_t = u_xx - F(u)_x * * * * *
% * * * * u(x,0)=sin(2*pi*x)*e^(-x^2)/2*delta
% * * * * Periodic boundary conditions u(a,t)=u(b,t)
% by METHOD OF LINES with centred second order finite difference scheme

clear all
close all

M=500;
a=0; b=6;
delta=0.5;
x=linspace(a,b,M)';
dx=(b-a)/(M-1); #space step discretization
A = toeplitz (sparse ([1, 1], [1, 2],[-2, 1] / dx^ 2, 1, M));
A(1,M)=1/dx^2;
A(M,1)=1/dx^2; #Periodic boundary conditions
b=ones(M,1);

%***** TIME INTEGRATION
T=.1;
t0=0;
ts=1000; #time steps
dt=(T-t0)/(ts-1);
u0=sin(2*pi*x).*exp((-x.^2)/(2*delta));
U=NaN(length(x),ts+1);
U(:,1)=u0;
u=u0;
t=t0;
for i=1:ts+1
  u=(speye(M)-dt*A)\(u-dt*b); #backward Euler
  U(:,i+1)=u; #update
  plot(x,u,'b-o')
  xlabel('x')
  ylabel('u(x,t*)')
  legend('u')
  t= t + dt;
  title(sprintf('Time= %0.3f',t));
  pause(0.00001);
end
