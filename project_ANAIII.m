  clear all
  close all

  %% Dati del problema

  m = 1000; %space steps
  x = linspace(0,6,m)';
  delta = .5;
  %u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum periodic
   u0 = exp(-(x-3).^2/(2*delta)); %initial datum neumann
  %u0 = 5*zeros(m,1);
%   u0 = 1.5*max(0,1-abs(x-3));
  v0 = ones(m,1);
  T = .1;
  dx = x(3)-x(2);
  %dt = dx./norm(x.^2,2); %pay attention to CFL condition
  
%% Caso 1
%  epsilon=realmin*ones(m,1);

%% Caso 2
% epsilon = (x<3)+0.1*(x>=3);

%% Caso 3
%  epsilon =0.5*ones(m,1);

%% Caso 4
%  epsilon = 1e-6*ones(m,1);

%% Caso 5
 epsilon =0.1 * ones(m,1);

dt=dx*10;

n = floor(T/dt)+1;

U = zeros(2*m,n);
U(:,1) = [u0;v0];

%% Implementazione del metodo

t =  0;
i = 1;

% %% Periodic Boundary Conditions
% figure;
% while t+dt<T
%   u = U(1:m,i);
%   v = U(m+1:end,i);
%   %Nuova v
%   %epsilon(end) = epsilon(1);
%   v(1:end-1) = epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*v(1:end-1)-...
%                 dt./(epsilon(1:end-1).^2+dt).*(([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)-(x(1:end-1).^2-6*x(1:end-1)).*u(1:end-1));
%   v(end) = v(1);
%   %Nuova u
%   u(1:end-1) = u(1:end-1) - dt*(epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*([v(2:end-1);v(1)]-[v(end-1);v(1:end-2)])/(2*dx) -...
%       dt./(epsilon(1:end-1).^2+dt).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 +...
%       dt./(epsilon(1:end-1).^2+dt).*(((x(1:end-1).^2-6*x(1:end-1)).*([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx))));
%   u(end) = u(1);
% 
%   
% % % Caso f(u)=u^2
% %   %Nuova v
% %   v(1:end-1) = epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*v(1:end-1)-dt./(epsilon(1:end-1).^2+dt).*(([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)-u(1:end-1).^2);
% %   v(end) = v(1);
% %   %Nuova u
% %   u(1:end-1) = u(1:end-1) - dt*(epsilon(1:end-1).^2./(epsilon(1:end-1).^2+dt).*([v(2:end-1);v(1)]-[v(end-1);v(1:end-2)])/(2*dx) -dt./(epsilon(1:end-1).^2+dt).*([u(2:end-1);u(1)]-2*u(1:end-1)+[u(end-1);u(1:end-2)])/dx^2 +dt./(epsilon(1:end-1).^2+dt).*(2*u(1:end-1).*([u(2:end-1);u(1)]-[u(end-1);u(1:end-2)])/(2*dx)));
% %   u(end) = u(1);
% 
%   %Aggiorno e memorizzo
%   U(1:m,i+1) = u;
%   U(m+1:end,i+1) = v;
%   i = i + 1;
%   t = t + dt;
%     
% end
% plot(x,u,'r-o',x,v,'b-o','MarkerSize',1)

%% Homogeneus Neumann Boundary Conditions

%Scrivendo u'(x1) = (u2-u0)/2h = 0, segue che u2 = u0, e analogamente
%un-1=un+1, quindi quando nello schema CFD c'è da usare u0 e un+1 mettiamo
%rispettivamente u2 e un-1, e lo risolviamo in tutti i nodi reali 1,..,n

figure;
while t+dt<T
  u = U(1:m,i);
  v = U(m+1:end,i);
  %Nuova v
  %epsilon(end) = epsilon(1);
  v(1:end) = epsilon(1:end).^2./(epsilon(1:end).^2+dt*ones(m,1)).*v(1:end)-...
                dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*(([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx)-(x(1:end).^2).*u(1:end));
  %Nuova u
  u(1:end) = u(1:end) - dt*(epsilon(1:end).^2./(epsilon(1:end).^2+dt*ones(m,1)).*([v(2:end);v(end-1)]-[v(2);v(1:end-1)])/(2*dx) -...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*([u(2:end);u(end-1)]-2*u(1:end)+[u(2);u(1:end-1)])/(dx^2) +...
      dt*ones(m,1)./(epsilon(1:end).^2+dt*ones(m,1)).*(((x(1:end).^2).*([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx))));

  
%   %Nuova v
%   %epsilon(end) = epsilon(1);
%   v(1:end) = epsilon(1:end).^2./(epsilon(1:end).^2+dt).*v(1:end)-...
%                 dt./(epsilon(1:end).^2+dt).*(([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx)-(u(1:end).^2));
%   %Nuova u
%   u(1:end) = u(1:end) - dt*(epsilon(1:end).^2./(epsilon(1:end).^2+dt).*([v(2:end);v(end-1)]-[v(2);v(1:end-1)])/(2*dx) -...
%       dt./(epsilon(1:end).^2+dt).*([u(2:end);u(end-1)]-2*u(1:end)+[u(2);u(1:end-1)])/(dx^2) +...
%       dt./(epsilon(1:end).^2+dt).*((2*u.*([u(2:end);u(end-1)]-[u(2);u(1:end-1)])/(2*dx))));
% 

  %Aggiorno e memorizzo
  U(1:m,i+1) = u;
  U(m+1:end,i+1) = v;
  i = i + 1;
  t = t + dt;
    
end
plot(x,u,'r-o',x,v,'b-o','MarkerSize',1)