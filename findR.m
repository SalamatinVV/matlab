%% 1
x= [0,1,2,3,4];
y= [0.03,0.98,2.09,2.98,4.11];
k=(sum(x.*y))/sum(x.^2);
x_real=[0:0.01:5];
y_real=x_real.*k;
plot(x,y, 'o');
hold on;
plot(x_real,y_real);
%% 1
x = [0,1,2,3,4];
y = [0.91,2.05,2.89,4.11,4.93];
sumSq= sum(x.*x);
sqSum=(sum(x))^2;
sumXY=sum(x.*y);
xySum=sum(y)*sum(x);
k=(sumXY-(1/length(y))*xySum)/(sumSq-(1/length(y))*sqSum);
b=(1/length(y))*sum(y)-(1/length(y))*sum(x)*k;
x_real=[0:0.01:5];
y_real=k*x_real+b;
plot(x,y, 'o');
hold on;
%% 2
x = [0,1,2,3,4];
y = [0.91,0.15,0.02,0.0025,0.0003];
g=0:0.1:5;
F=zeros(1,length(g));
A=zeros(1,length(g));
for k=1:length(g)
    C1=sum(exp((-g(k))*x).*y);
    C2=sum(exp((-2*g(k))*x));
    C3=sum(exp((-2*g(k))*x).*x);
    C4=sum(exp((-g(k))*x).*(x.*y));
    F(k)=C1*C3-C2*C4;
    for i=2:length(g)
        if F(i)*F(i-1)<0
            rootgamma=(g(i)+g(i-1))/2;
        end
    end
    A=sum((exp((-rootgamma)*x).*y)/sum(exp((-2*rootgamma)*x)));
    x_real=[0:0.01:5];
    y_real=A*exp(-rootgamma.*x_real);
end
plot(x,y,'o');
hold on;
plot(x_real,y_real);
%% 3
x = [0,1,2,3,4,5,6];
y = [1.5,1.8,2,2.1,1.4,0.9,0.7];
B=y.';
N=7;
A=zeros(N);
for i=1:N
    for j=1:N
        A(i,j)=(x(i))^(j-1);
        
    end
end
a=A^(-1)*B;
plot(x,y);
%% 3 square spline 
x = 0:1:100;
y = rand(1, length(x));
c = zeros(1,length(x));
c(1) = 0;
h=1;
for i=2:length(y)-1
    c(i)=(1/(h^2))*(y(i+1)-2*y(i)+y(i-1))-c(i-1);
end
xarr = [];
yarr = [];
for j=1:length(y)-1
    xf=x(j):(x(j+1)-x(j))/100:x(j+1);
    yf=y(j)+(y(j+1)-y(j))*(xf-x(j))/(x(j+1)-x(j))+c(j)*(xf-x(j)).*(xf-x(j+1));
    xarr = [xarr xf];
    yarr = [yarr yf];
end
plot(x,y,'o');
hold on;
plot (xarr,yarr);

%% 4 cub spline
x = [ 0 1 2 3 4 5 6];
y = [1.5 2 2.2 1.8 1.5 0.9 0.7];
c = zeros(1,length(x));
d = zeros(1,length(x));
c(2) = 0;
d(2)=0;
h=1;
for i=2:length(y)-1
    d(i+1)=d(i)-6*c(i)/h +6*(y(i+1)-2*y(i)+y(i-1))/h^3;
    c(i+1)=c(i)+h*d(i+1);
end
xarr = [];
yarr = [];
for j=2:length(y)-1
    xf=x(j):(x(j+1)-x(j))/7:x(j+1);
    yf=y(j)+(((y(j)-y(j-1))/h)+h*c(j)/2-h^2*d(j)/6)*(xf-x(j))+(((xf-x(j)).^2)*c(j)/2)+(((xf-x(j)).^3)*d(j))/6;
    xarr = [xarr xf];
    yarr = [yarr yf];
end

plot(x,y,'o');
hold on;
plot (xarr,yarr);

%% 5 ÑËÀÓ 
A=[0 0 3 2 7; 0 1 1 3 5; 6 8 0 1 7; 3 3 2 8 5; 9 6 7 0 1];
b=[1;5;7;3;2];
n=length(b); 
reffA=rref(A);
for k=1:n-1
    if A(k,k)==0
        for i=k+1:n
            if A(i,k)~=0
                temp=A(k,:);
                A(k,:)=A(i,:);
                A(i,:)=temp;
                temp=b(k);
                b(k)=b(i);
                b(i)=temp;
                break
            end
        end
    end
end
if sum((A(:,1)).^2)~=0
    for i=2:n
        b(i)=b(i)+(-A(i,1)/A(1,1))*b(1);
        A(i,:)=A(i,:)+(-A(i,1)/A(1,1)*A(1,:));
    end
end
if sum((A(2:end,2)).^2)~=0
    for i=2:n
        if A (i,2)~=0
            rem = A(i,:);
            remnum=i;
            remb=b(i);
        end
    end
    A(remnum,:)=A(2,:);
    A(2,:)=rem;
    b(remnum)=b(2);
    b(2)=remb;
    for i=3:n
        b(i)=b(i)+(-A(i,2)/A(2,2))*b(2);
        A(i,:)=A(i,:)+(-A(i,2)/A(2,2)*A(2,:));
    end
end

%% 6 

b= [7;5;3;2];
a=[0;2;0;5];
c=[9;8;1;0];
f=[7;1;4;2];
N=length(f);
beta=zeros(N-1,1);
alpha=zeros(N-1,1);
alpha(1)=-c(1)/b(1);
beta(1)=f(1)/b(1);
for i=2:N-1
    alpha(i)=-c(i)/(a(i)*alpha(i-1)+b(i));
    beta(i)=(f(i)-a(i)*beta(i-1))/(a(i)*alpha(i-1)+b(i));
end
x=zeros(N,1);
x(N)=(f(N)-a(N)*beta(N-1))/(a(N)*alpha(N-1)+b(N));
for i=N-1:-1:1
    x(i)=alpha(i)*x(i+1)+beta(i);
end

%% 7
lim1=0.01;
lim2=10;
N=100;
x=lim1:(lim2-lim1)/N:lim2;
y=zeros(length(x));
y=sin(x)./x-0.5;
M=10;
for k = 1:M
    x=lim1:(lim2-lim1)/N:lim2;
    y=sin(x)./x-0.5;
    for i=1:length(x)-1
        if y(i)*y(i+1)<0
          lim1=x(i);
          lim2=x(i+1);
        end
    end
end
root=(lim1+lim2)/2;

%%

lim1=[2.5;6.58;8.3];
lim2=[3;7.2;8.5];
x=zeros(3,101);
N=100;
y=zeros(length(x));
y=sin(x)./x-0.1;
y1=zeros(length(x));
M=10;
for k = 1:M
    for k = 1:length(lim1)
        x(k,:)=lim1(k):(lim2(k)-lim1(k))/N:lim2(k);
    end
    y=sin(x)./x-0.1;
    for k =1:length(lim1)
    for i=1:length(x)-1
        if y(k,i)*y(k,i+1)<0
          lim1(k)=x(k,i);
          lim2(k)=x(k,i+1);
        end
    end
    end
end
root=(lim1+lim2)/2;







    