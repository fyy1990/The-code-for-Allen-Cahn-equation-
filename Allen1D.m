% Allen-Cahn equation
%u_{t}=\epsion^2\laplace*u+u-u^3;
%(I-frac{\tau \epsion^2}{2}D_{h}-frac{\tau}{2})U^{n+1}=(I+frac{\tau
%\epsion^2}{2}D_{h}+frac{\tau}{2})U^{n}-((U^{n})^2+(U^{n+1})^2)*(U^n+U^{n+1})/4;数值格式；
clear;
%% Time and space grid partitioning
h=1/4;% space step
alpha=1.2;
x=0:h:1;x=x';
n=length(x);
x1=x(2:end-1);
u0=x1.^5.*(1-x1).^5;% initial value
epsion=0.01;
dt=0.001;%time step
tend=1;
T=0:dt:tend;
N=length(T);
%% difference matrix
b=zeros(1,n-2);
for i=1:n-2
    b(i)=(-1)^(i-1)*gamma(alpha+1)/(gamma(alpha/2-(i-1)+1)*gamma(alpha/2+i-1+1));
end

g=zeros(1,n-2);
for kk=1:n-2
if mod(kk+1, 2)== 0
    
   g(kk)=4/3*b(kk)-1/(3*2^(alpha))*b((kk+1)/2);% number is even
else
   g(kk)=4/3*b(kk);% number is odd
end
end
cc=fliplr(g(2:end))';
a=[g';0;cc];
diagg=fft(a);
%% Right end term f
  for i=1:N
  f(:,i)=exp(-3*(T(i)+dt/2)).*x1.^15.*(1-x1).^15-2*exp(-(T(i)+dt/2)).*x1.^5.*(1-x1).^5+...
        epsion^2/(2*cos(alpha*pi/2))*exp(-(T(i)+dt/2)).*(gamma(6)/gamma(6-alpha)*((x1).^(5-alpha)+(1-x1).^(5-alpha))-...
    5*gamma(7)/gamma(7-alpha)*((x1).^(6-alpha)+(1-x1).^(6-alpha))+...
            10*gamma(8)/gamma(8-alpha)*((x1).^(7-alpha)+(1-x1).^(7-alpha))-10*gamma(9)/gamma(9-alpha)*(x1.^(8-alpha)+(1-x1).^(8-alpha))+...
        5*gamma(10)/gamma(10-alpha)*((x1).^(9-alpha)+(1-(x1)).^(9-alpha))-gamma(11)/gamma(11-alpha)*((x1).^(10-alpha)+(1-(x1)).^(10-alpha)));      
     
end
   uex=exp(-T(end))*(x1).^5.*(1-x1).^5;
%% Nonlinear iterative format calculation
U=zeros(N-1,1);yy=zeros(n-2,1); u1=zeros(n-2,1);  
U(1)=max(abs(u0));
for j=1:N-1
   tol = 1; uu = u0; count = 0;dd=dt/2*(f(:,j));
   a1=dt/8;b1=dt/8*u0;c1=dt/8*u0.^2+1/2*(3/2-dt/2);
    while tol > 1e-10
        
        rhss=dt*epsion^2/(2*h^(alpha))*(ifft(diagg.*fft([(u0+uu)/2;zeros(n-2,1)])));
        rhs=1/2*(u0+uu)/2-rhss(1:n-2)+u0+dd;
        
        d1=(3/2-dt/2)/2*u0+dt/8*u0.^3-rhs;

        for k=1:n-2
           xx=roots([a1 b1(k) c1(k) d1(k)]);
            for l=1:3
                xxx(l)=isreal(xx(l));
            end
            u1(k)=xx(xxx);
        end             
    
       tol = max(abs(u1-uu));
        uu = u1; 
    end
    
    u0 = u1; 
    
    U(j+1)=max(abs(u0));
    
   
     %t = j*dt;
  
%     plot(x,abs(u0)), drawnow
    
  %  fprintf('Computing times = %d    Total times = %d\n',j,N)
     
end

 err=max(abs(u0-uex))

C=ones(N,1);
figure(1)
plot(T,U,'b--o','LineWidth', 2.5)
hold on
plot(T,C,'--r','LineWidth', 2.5)
hl=legend('\epsilon=0.1,\tau=1.5','C=1');
set(hl,'Box','off','Location','east','FontSize',15);
xlabel('t')
ylabel('|u|_{max}')


































