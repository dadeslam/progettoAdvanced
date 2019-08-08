%%% Solution to | u_t + v_x=0
%               |
%               |v_t + 1/eps^2 u_x=-1/eps^2(v-u)
%               |
%               |u(x,0)=                    x \in [-L,L], t>0
%               |v(x,0)=
%               | + Periodic Boundary Conditions

%Analytical solutions: |u(x,t)=exp(-try)*sin(x-t)
%                      |v(x,t)=exp(-t)*(sin(x-t)-cos(x-t))
%We choose T=1, epsilon=1e-6, and discretize space derivatives with second order finite differences

clear all
close all
L=6;
m = 3000; %space steps
x = linspace(0,L,m)';
T = .5; %final time
dx = x(3)-x(2);

%u0, v0 the ones of the paper
% u0=sin(x);
% v0=sin(x)-cos(x);
delta=.5;
u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum
% v0=zeros(m,1);
u0=cos(pi*x); %Neumann boundary condition
v0=zeros(m,1);

fprintf('Case 1: realmin \n')
fprintf('Case 2: 0.1 \n' )
fprintf('Case 3: 0.5 \n' )
fprintf('Case 4: 1 \n' )
fprintf('Case 5: variable in symmetric\n' )
fprintf('Case 6: variable \n' )
n=input('Enter a case for epsilon: ');

switch n
  case 1
    epsilon=realmin*ones(m,1);
  case 2
    epsilon =0.1*ones(m,1);
  case 3
    epsilon =0.5*ones(m,1);
  case 4
    epsilon =1*ones(m,1);
  case 5
    epsilon = (x<0)+0.1*(x>=0);
  case 6
    epsilon = (x<3)+0.1*(x>=3);
end

dt=dx/50; %CFL condition
n = floor(T/dt)+1;
t =  0;
i = 1;

%% Periodic Boundary Conditions
figure('units','normalized','outerposition',[0 0 1 1]) %run immediately full screen
u=u0;
v=v0;
while t+dt<T %main loop
%   %New v
%   vold=v;
%   v(1:end-1) = epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*v(1:end-1)-...
%                 dt./(epsilon(1:end-1).^2+dt).*(([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)-(x(1:end-1).^2).*u(1:end-1));
%   v(end) = v(1);
%   %New u
%   u(1:end-1) = u(1:end-1) - dt*(epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*([vold(2:end-1);vold(1)]-[vold(end-1);vold(1:end-2)])/(2*dx) -...
%       dt./(epsilon(1:end-1).^2+dt).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 +...
%       dt./(epsilon(1:end-1).^2+dt).*(((x(1:end-1).^2).*([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx))));
%   u(end) = u(1);


%   %New v

   %Stack suggestion
  vold=v;
  v(1:end) = epsilon(1:end).^2./(epsilon(1:end).^2+dt*ones(m,1)).*v(1:end)-...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*(([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx)-(x(1:end).^2).*u(1:end));
  v(1)=0;
  v(end)=0;
  %Nuova u
  u(1:end) = u(1:end) - dt*(epsilon(1:end).^2./(epsilon(1:end).^2+dt*ones(m,1)).*([2*vold(2);vold(3:end);2*vold(end)]-[2*vold(1);vold(1:end-2);2*vold(end-1)])/(2*dx) -...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*([u(2:end);u(end-1)]-2*u(1:end)+[u(2);u(1:end-1)])/(dx^2) +...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*(((x(1:end).^2).*([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx))));
  


% plot(x,u,'b',x,v,'r')
% title(sprintf('Comparison at Time= %0.3f',t));
% legend('u','v');
% grid on
% drawnow

  i = i + 1;
  t = t + dt;
end
plot(x,u,'b',x,v,'r')
legend('u','v');
%axis([4,6,-2,2])
grid on
