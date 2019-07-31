  clear all
  close all

  %% Dati del problema

  m = 400; %space steps
  x = linspace(0,6,m)';
  T = .1;
  dx = x(3)-x(2);
  dt = dx/m; %pay attention to CFL condition
  n = floor(T/dt)+1;
%   %% Caso 1
%   vec=realmin;



% %% Caso 2
% vec = [1;0.5;0.1];
% vec = vec(1);

delta = .5;
u0 = sin(2*pi*x).*exp(-x.^2/(2*delta)); %initial datum

v0 = zeros(m,1);

F = @(t) t;

U = zeros(2*m,n);
U(:,1) = [u0;v0];

%% Implementazione del metodo

for epsilon = vec
    t = 0;
    i = 1;
    %figure;
    while t+dt<T
      u = U(1:m,i);
      v = U(m+1:end,i);
      %Nuova v
      v(1:end) = epsilon^2/(epsilon^2+dt)*v(1:end) - dt/(epsilon^2+dt)*( ([u(2:end);u(2)]-[u(end-1);u(1:end-1)])/(2*dx) - F(u(1:end)) );
      %Nuova u
      u(1:end) = u(1:end) - dt*(epsilon^2/(epsilon^2+dt)*([v(2:end);v(2)]-[v(end-1);v(1:end-1)])/(2*dx) -...
          dt/(epsilon^2+dt)*([u(2:end);u(2)]-2*u(1:end)+[u(end-1);u(1:end-1)])/dx^2 + ...
          dt/(epsilon^2+dt)*([u(2:end);u(2)]-[u(end-1);u(1:end-1)])/(2*dx) );
      %Aggiorno e memorizzo
      U(1:m,i+1) = u;
      U(m+1:end,i+1) = v;
      i = i + 1;
      t = t + dt;
      %plot(x,u,'r--',x,v,'b--','MarkerSize',1)
      plot(x,u,'r--');
      axis([0 6 -.5 .5])
      xlabel('x')
      ylabel('u(x,t*)')
      %legend('u','v')
      title(sprintf('Time= %0.3f',t));
      pause(0.0000000001);
    end
end


