clear all

%Example 1D
% delta=0.01;
% del=-2:delta:-0.9;
% dell(1)=-2;
% N=length(del);
% y(1)=(-1-del(1))^(1/2);

%Example 2D====
delta=0.001;
del=-1:delta:-0.;
dell(1)=-1;
N=length(del);
y(1)=(-1-2*del(1))^(1/2);

%%===For function y^2+1+\delta=0===
%Eg 1D
%f=@(x,y)-1/(2*y);

%Eg 2D
f=@(x,y)-1/y;
g=@(x,y)(-1-x-y)^(1/2);
%%
%%===Euler Method====
for i=1:N
    y(i+1)=y(i)+f(del(i),y(i));
    sol=g(del(i),del(i));
end
figure(4)
plot(del,y)
hold on
plot(del,sol)


%%
%%===Modified Euler method=====

for i=1:N-1
    y(i+1)=y(i)+delta*f(del(i),y(i));
    y(i+1)=y(i)+1/2*delta*(f(del(i),y(i))+f(del(i+1),y(i+1)));
    
    sol(i)=g(del(i),del(i));
end

figure(1)
plot(del,y)
hold on
plot(del(1:N-1), sol,'r--')
legend('Approx','Exact')
%plot(del(1:ind),(-1-del(1:ind)).^(1/2))
%%
%%======Third Order Runge Kutta Method====
for i=1:N
    k1 = -1/(2*y(i));
    k2 = -1/(2*(y(i)+1/2*delta*k1));
    k3 = -1/2*1/(y(i)-delta*k1+2*delta*k2);
    y(i+1) = y(i)+delta/6*(k1+4*k2+k3);
    dell(i+1)=dell(i)+delta;
end
ind=find(imag(y)>0, 1 );
figure(2)
plot(dell,y)
hold on
plot(del(1:ind),(-1-del(1:ind)).^(1/2))

%%
%%====Classic Fourth Order R-K Method====
for i =1:N
    k1=f(del(i),y(i));
    k2=f(del(i)+1/2*delta,y(i)+delta*1/2*k1);
    k3=f(del(i)+1/2*delta,y(i)+delta*1/2*k2);
    k4=f(del(i)+delta,y(i)+delta*k3);
    y(i+1)=y(i)+delta*(k1/6+k2/3+k3/3+k4/7);
        dell(i+1)=dell(i)+delta;
end
ind=find(del==-1);
figure(3)
plot(dell,y)
hold on
plot(del(1:ind),(-1-del(1:ind)).^(1/2))

%%
%%==Puiseux serious with three terms====
a0=-8.580022165053089e-15; 
a1=1.000000000000009e+00 ;
a2=-2.704753515664180e-15 ;
c=-9.999999999999947e-01 ;
t0=1.732050807568879e+00 ;
t1=1.414213562373097e+00 ;
t2=1.000000000000003e+00 ;
t3=7.071067811865516e-01 ;
d0=-4;
d1=-3;d2=-2;d3=-1.5;
clear y
for i=1:N
    y(i)=a0+a1*(c-del(i))^(1/2)+a2*(c-del(i));
end
figure(5)
plot(del,y)
hold on
plot(del,(-1-del).^(1/2))
    