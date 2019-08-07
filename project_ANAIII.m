  %% Dati del problema

  m = 400; %space steps
  m = 1000; %space steps
  x = linspace(0,6,m)';
  delta = .5;
  %u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum
  % u0 = 5*ones(m,1);
  u0 = 1.5*max(0,1-abs(x-3));
  v0 = zeros(m,1);
  %u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum periodic
   u0 = exp(-(x-3).^2/(2*delta)); %initial datum neumann
  %u0 = 5*zeros(m,1);
%   u0 = 1.5*max(0,1-abs(x-3));
  v0 = ones(m,1);
  T = .1;
  dx = x(3)-x(2);
  %dt = dx./norm(x.^2,2); %pay attention to CFL condition

%% Caso 1
% epsilon=realmin*ones(m,1);
%  epsilon=realmin*ones(m,1);

%% Caso 2
% epsilon = (x<3)+0.1*(x>=3);

%% Caso 3
%epsilon =ones(m,1);
%  epsilon =0.5*ones(m,1);

%% Caso 4
 epsilon = 1e-6*ones(m,1);
%  epsilon = 1e-6*ones(m,1);

%% Caso 5
% epsilon = (x<3-dx)+0.1*(x>=3 & x<=6-dx)+(-0.9/dx*(x-3+dx)+1).*(x>=3-dx & x<3)+(0.9/dx*(x-6+dx)+0.1).*(x>6-dx);
% epsilon(end) = epsilon(1);

dt=dx/100;
 epsilon =0.1 * ones(m,1);

dt=dx*10;

n = floor(T/dt)+1;

@@ -39,8 +39,8 @@

%% Implementazione del metodo

% t =  0;
% i = 1;
t =  0;
i = 1;

% %% Periodic Boundary Conditions
% figure;
@@ -76,46 +76,41 @@
% end
% plot(x,u,'r-o',x,v,'b-o','MarkerSize',1)

%% Homogeneus Neumann Boundary Conditions

%Scrivendo u'(x1) = (u2-u0)/2h = 0, segue che u2 = u0, e analogamente
%un-1=un+1, quindi quando nello schema CFD c'è da usare u0 e un+1 mettiamo
%rispettivamente u2 e un-1, e lo risolviamo in tutti i nodi reali 1,..,n

%% Outflow Boundary Conditions
t=0;
i=1;
figure;
while t+dt<T
  u = U(1:m,i);
  v = U(m+1:end,i);
  %Nuova v
  %epsilon(end) = epsilon(1);
  v(2:end-1) = epsilon(2:end-1).^2./(epsilon(2:end-1).^2+dt).*v(2:end-1)-...
                dt./(epsilon(2:end-1).^2+dt).*((u(3:end)-u(1:end-2))/(2*dx)-u(2:end-1));
  v(1) = v(2);
  v(end) = v(end-1);

  v(1:end) = epsilon(1:end).^2./(epsilon(1:end).^2+dt*ones(m,1)).*v(1:end)-...
                dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*(([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx)-(x(1:end).^2).*u(1:end));
  %Nuova u
  u(2:end-1) = u(2:end-1) - dt*(epsilon(2:end-1).^2./(epsilon(2:end-1).^2+dt).*(v(3:end)-v(1:end-2))/(2*dx) -...
      dt./(epsilon(2:end-1).^2+dt).*(u(3:end)-2*u(2:end-1)+u(1:end-2))/dx^2 +...
      dt./(epsilon(2:end-1).^2+dt).*((u(3:end)-u(1:end-2))/(2*dx)));
  u(1) = u(2);
  u(end) = u(end-1);

  %Da riadattare alle nuove cond bordo
% % Caso f(u)=u^2
  u(1:end) = u(1:end) - dt*(epsilon(1:end).^2./(epsilon(1:end).^2+dt*ones(m,1)).*([v(2:end);v(end-1)]-[v(2);v(1:end-1)])/(2*dx) -...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*([u(2:end);u(end-1)]-2*u(1:end)+[u(2);u(1:end-1)])/(dx^2) +...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*(((x(1:end).^2).*([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx))));


%   %Nuova v
%   v(1:end-1) = epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*v(1:end-1)-dt./(epsilon(1:end-1).^2+dt).*(([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)-u(1:end-1).^2);
%   v(end) = v(1);
%   %epsilon(end) = epsilon(1);
%   v(1:end) = epsilon(1:end).^2./(epsilon(1:end).^2+dt).*v(1:end)-...
%                 dt./(epsilon(1:end).^2+dt).*(([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx)-(u(1:end).^2));
%   %Nuova u
%   u(1:end-1) = u(1:end-1) - dt*(epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*([v(2:end-1);v(1)]-[v(end-1);v(1:end-2)])/(2*dx) -dt./(epsilon(1:end-1).^2+dt).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 +dt./(epsilon(1:end-1).^2+dt).*(2*u(1:end-1).*([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)));
%   u(end) = u(1);
%   u(1:end) = u(1:end) - dt*(epsilon(1:end).^2./(epsilon(1:end).^2+dt).*([v(2:end);v(end-1)]-[v(2);v(1:end-1)])/(2*dx) -...
%       dt./(epsilon(1:end).^2+dt).*([u(2:end);u(end-1)]-2*u(1:end)+[u(2);u(1:end-1)])/(dx^2) +...
%       dt./(epsilon(1:end).^2+dt).*((2*u.*([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx))));
% 

  %Aggiorno e memorizzo
  U(1:m,i+1) = u;
  U(m+1:end,i+1) = v;
  i = i + 1;
  t = t + dt;
  %plot(x,u,'r-o',x,v,'b-o','MarkerSize',1)
  %pause(0.00001);
%   xlabel('x')
%   ylabel('u(x,t*)')
%   legend('u','v')
%   title(sprintf('Time= %0.3f',t));

end
plot(x,u,'r-o',x,v,'b-o','MarkerSize',1) 
