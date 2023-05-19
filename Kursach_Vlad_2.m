h=0.01;
x=1:h:2;
N=length(x);
y=zeros(1,N);
a=zeros(1,N);
b=zeros(1,N);
c=zeros(1,N);
f=zeros(1,N);

a(1)=0;
b(1)=1-1/h;
c(1)=1/h;
f(1)=0;
for i=2:N-1
a(i)=1/(h^2);
b(i)=-2/h^2+(x(i))^3*sin(x(i))/h+(x(i)^2)*exp(-2*x(i));
c(i)=1/(h^2)-x(i)^3*sin(x(i))/h;
f(i)=sqrt(log10(x(i)+2));
end
a(N)=2/h;
b(N)=1 - 2/h;
f(N)=4;

alpha=zeros(1,N);
beta=zeros(1,N);
alpha(1)=-c(1)/b(1);
beta(1)=f(1)/b(1);

for i=2:N-1
    alpha(i)=-c(i)/(a(i)*alpha(i-1)+b(i));
    beta(i)=(f(i)-a(i)*beta(i-1))/(a(i)*alpha(i-1)+b(i));
end
y(N)=(f(N)-a(N)*beta(N-1))/(a(N)*alpha(N-1)+b(N));
for i=N-1:-1:1
    y(i)=alpha(i)*y(i+1)+beta(i);
end
plot(x,y);

