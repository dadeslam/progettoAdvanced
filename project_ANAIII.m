%%% Solution to | u_t + v_x=0
%               |
%               |v_t + 1/eps^2 u_x=-1/eps^2(v-u)
%               |
%               |u(x,0)=sin(x)                     x \in [-pi,pi], t>0
%               |v(x,0)=sin(x)-cos(x)
%               | + Periodic Boundary Conditions

%Analytical solutions: |u(x,t)=exp(-try)*sin(x-t)
%                      |v(x,t)=exp(-t)*(sin(x-t)-cos(x-t))
%We choose T=1, epsilon=1e-6, and discretize space derivatives with second order finite differences

clear all
close all
m = 100; %space steps
x = linspace(-pi,pi,m)';
T = 1; %final time
dx = x(3)-x(2);

%u0, v0 the ones of the paper
u0=sin(x);
v0=sin(x)-cos(x);
%% Caso 1
%   epsilon=realmin*ones(m,1);
%% Caso 2
%epsilon = (x<3)+0.1*(x>=3);
%% Caso 3
%epsilon =0.1*ones(m,1);
%% Caso 4
epsilon = 1e-6*ones(m,1); %the one of the paper
%% Caso 5
%  epsilon =0.1 * ones(m,1);

dt=dx/50; %CFL condition
n = floor(T/dt)+1;
t =  0;
i = 1;

%% Periodic Boundary Conditions
figure('units','normalized','outerposition',[0 0 1 1]) %run immediately full screen
u=u0;
v=v0;
while t+dt<T %main loop
  %New v
  vold=v;
  v(1:end-1) = epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*v(1:end-1)-...
                dt./(epsilon(1:end-1).^2+dt).*(([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)-(x(1:end-1).^0).*u(1:end-1));
  v(end) = v(1);
  %New u
  u(1:end-1) = u(1:end-1) - dt*(epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*([vold(2:end-1);vold(1)]-[vold(end-1);vold(1:end-2)])/(2*dx) -...
      dt./(epsilon(1:end-1).^2+dt).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 +...
      dt./(epsilon(1:end-1).^2+dt).*(((x(1:end-1).^0).*([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx))));
  u(end) = u(1);


plot(x,u,'b-d',x,v,'k-d',x,exp(-t)*sin(x-t),'r-o',x,exp(-t)*(sin(x-t)-cos(x-t)),'y-o')
title(sprintf('Comparison at Time= %0.3f',t));
legend('u','v','u_{ex}','v_{ex}');
grid on
drawnow

  i = i + 1;
  t = t + dt;
end
