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
m = 201; %space steps
x = linspace(0,L,m)';
T = 3; %final time
dx = x(3)-x(2);

%u0, v0 the ones of the paper
% u0=sin(x);
% v0=sin(x)-cos(x);
delta=.5;
%u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum
% v0=zeros(m,1);
%u0=cos(pi*x); %Neumann boundary condition
u0 = sin(4*pi*(x-3)).*exp(-4*(x-3).^2/(2*delta)); %WORKS
%u0 = x.^2.*(x-1).*(x-6).^2/400; %WORKS
%u0=ones(m,1);
v0=zeros(m,1);

fprintf('Case 1: 1e-6 \n')
fprintf('Case 2: 0.1 \n' )
fprintf('Case 3: 0.5 \n' )
fprintf('Case 4: 1 \n' )
fprintf('Case 5: variable \n' )
n=input('Enter a case for epsilon: ');

pos = 1;
switch n
  case 1
    epsilon=1e-6.*ones(m,1);
  case 2
    epsilon =0.1*ones(m,1);
  case 3
    epsilon =0.5*ones(m,1);
  case 4
    epsilon =1*ones(m,1);
  case 5
    epsilon = (x<=3)+0.1*(x>3);
    pos = find(diff(epsilon));
    
end

dt=dx^2/400; %CFL condition
n = floor(T/dt)+1;
t =  0;
i = 1;

%% Periodic Boundary Conditions
figure('units','normalized','outerposition',[0 0 1 1]) %run immediately full screen
u=u0;
v=v0;
while t+dt<T %main loop
    %% Funziona con epsilon costante, dt=dx^2/400 e m = 201 
    if(pos==1) %Epsilon is constant in this case
        %New v
      vold=v;
      v(1:end-1) = epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*v(1:end-1)-...
                    dt./(epsilon(1:end-1).^2+dt).*(([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)-((x(1:end-1).^2-6*x(1:end-1))).*u(1:end-1));

      %New u
      u(1:end-1) = u(1:end-1) - dt*(epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*([vold(2:end-1);vold(1)]-[vold(end-1);vold(1:end-2)])/(2*dx) -...
          dt./(epsilon(1:end-1).^2+dt).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 +...
          dt./(epsilon(1:end-1).^2+dt).*(((x(1:end-1).^2-6*x(1:end-1)).*([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx))));
      
      v(end) = v(1);
      u(end) = u(1);
    
    else %Here the epsilon is not constant
       %New v
      vold=v;
      
      v(1:pos) = epsilon(1:pos).^2./(epsilon(1:pos).^2+dt).*v(1:pos)-...
                    dt./(epsilon(1:pos).^2+dt).*(([u(2:pos);u(pos+1)]-[u(end-1);u(1:pos-1)])/(2*dx)-(x(1:pos).^2-6*x(1:pos)).*u(1:pos));
      u(1:pos) = u(1:pos) - dt*(epsilon(1:pos).^2./(epsilon(1:pos).^2+dt).*([vold(2:pos);vold(pos+1)]-[vold(end-1);vold(1:pos-1)])/(2*dx) -...
          dt./(epsilon(1:pos).^2+dt).*([u(2:pos);u(pos+1)]-2*u(1:pos)+[u(end-1);u(1:pos-1)])/dx^2 +...
          dt./(epsilon(1:pos).^2+dt).*((x(1:pos).^2 -6*x(1:pos)).*([u(2:pos);u(pos+1)]-[u(end-1);u(1:pos-1)])/(2*dx))); 
            
      v(pos+1:end-1) = epsilon(pos+1:end-1).^2./(epsilon(pos+1:end-1).^2+dt).*v(pos+1:end-1)-...
                    dt./(epsilon(pos+1:end-1).^2+dt).*((u(pos+2:end)-u(pos:end-2))/(2*dx)-(x(pos+1:end-1).^2-6*x(pos+1:end-1)).*u(pos+1:end-1));
      u(pos+1:end-1) = u(pos+1:end-1) - dt*(epsilon(pos+1:end-1).^2./(epsilon(pos+1:end-1).^2+dt).*(vold(pos+2:end)-vold(pos:end-2))/(2*dx) -...
          dt./(epsilon(pos+1:end-1).^2+dt).*(u(pos+2:end)-2*u(pos+1:end-1)+u(pos:end-2))/dx^2 +...
          dt./(epsilon(pos+1:end-1).^2+dt).*(((x(pos+1:end-1).^2-6*x(pos+1:end-1)).*(u(pos+2:end)-u(pos:end-2))/(2*dx))));
      
      v(end) = v(1);
      u(end) = u(1);
    end
   i = i + 1;
  t = t + dt;
end
plot(x,u,'b',x,v,'r')
hold on;
plot(x(pos),u(pos),'go',x(pos),v(pos),'kx')
legend('u','v');
grid on