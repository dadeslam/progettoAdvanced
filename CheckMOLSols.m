% * * * * |u_t + v_x= 0 * * * * *
% * * * * |v_t + 1/eps^2 u_x = -1/eps^2(v-F(u))
% * * * * |u(x,0)=sin(4*pi*(x-3)).*exp(-4*(x-3).^2)
% * * * * |v(x,0)=0
% * * * * Periodic boundary conditions

% by METHOD OF LINES with centred second order finite difference scheme

clear all
close all

m=1000;
a=0; b=6;
delta=0.5;
x=linspace(a,b,m)';
dx=(b-a)/(m-1); %space step discretization

A1= toeplitz(sparse(1,2,-1/(2*dx),1,m), sparse(1,2,1/(2*dx),1,m));
A2= toeplitz(sparse(1,2,-1/(2*dx),1,m), sparse(1,2,1/(2*dx),1,m));


epsilon=@(x) 1*(x<=3) + 0.1*(x>3);
%epsilon=@(x) ones(m,1);
EPS=repmat(epsilon(x).^2,1,m); %epsilon^2 matrix
A2=A2./EPS;

%Build the system
u0=sin(4*pi*(x-3)).*exp(-4*(x-3).^2); 
%u0=x.^2.*(x-1).*(x-6).^2/400;
v0=zeros(m,1);

Z=zeros(m,m);
B=[Z,A1;A2,Z];
%f=@(y) (x.^2 - 6*x).*u;
f=@(y) (x.^2 - 6*x).*y(1:m);
%F=@(y) [zeros(m,1);0;(-1./epsilon(x(2:m-1)).^2).*(y(m+1:2*m)-f(y(1:m)))(2:m-1);0]; %boundary Conditions in Octave

V=[0;ones(m-2,1);0];
F=@(y) [zeros(m,1);((-1./epsilon(x(1:m)).^2).*(y(m+1:2*m)-f(y(1:m)))).*V]; %boundary Conditions
%Boundary Conditions

%periodic in v
A1(1,1:2)=[0,1]/dx;
A1(1,m)=-1/dx;

A1(m,1)=1/dx;
A1(m,m-1:m)=[-1,0]/dx;


%periodic in u
A2(1,1:2)=[0,1]/dx;
A2(1,m)=-1/dx;

A2(m,1)=1/dx;
A2(m,m-1:m)=[-1,0]/dx;

%% Solution of the ODEs system

y0=[u0;v0];
tspan=[0,1];
odefun=@(t,y) -B*y +F(y);


[t,y]=ode45(odefun,tspan,y0);
figure('units','normalized','outerposition',[0 0 1 1]) %run immediately full screen

for j=1:length(t)
  plot(x,y(j,1:m),'b',x,y(j,m+1:2*m),'r')
  xlabel('x')
  ylabel('')
  title(sprintf('Time= %0.3f',t(j)));
  drawnow
end
