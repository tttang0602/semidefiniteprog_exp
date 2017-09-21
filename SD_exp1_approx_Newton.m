% d0=-0.489;
% d1=-0.49;
% d2=-0.491;
% d3=-0.492;
% d4=-0.493;
% y0=0.204305518472662;
% y1=0.201107684592419;
% y2=0.197482095669551;
% y3=0.193115453265627;
% y4=0.186900321003191;% 

% d0=-0.49;
% d1=-0.39;
% d2=-0.29;
% d3=-0.19;
% d4=-0.09;
% 
% y0=0.201107684592418;
% y1=0.361285709535934;
% y2=0.499692256984107;
% y3=0.640735499050238;
% y4=0.785063583433010;

%%Four data points for Ex2 with fixed delta1 and delta2 starting at -0.05
%%and step size 0.005
%load('Ex2_sol_fine.mat')
load('Ex2_sol.mat')
% d0=-0.05;
% d1=-0.049;
% d2=-0.048;
% d3=-0.047;
% d4=-0.046;
% y0=sol_0_del2(1,1,1);
% y1=sol_0_del2(1,2,1);
% y2=sol_0_del2(1,3,1);
% y3=sol_0_del2(1,4,1);
% y4=sol_0_del2(1,5,1);
% 

d0=0.03;
d1=0.035;
d2=0.04;
d3=0.045;
d4=0.05;

% d0=-0.05;
% d1=-0.045;
% d2=-0.04;
% d3=-0.035;
% d4=-0.03;
ind2=11;
y0=sol(17,ind2,1);
y1=sol(18,ind2,1);
y2=sol(19,ind2,1);
y3=sol(20,ind2,1);
y4=sol(21,ind2,1);
% d0=-0.05;
% d1=-0.045;
% d2=-0.04;
% d3=-0.035;
% d4=-0.03;
% y0=sol(1,1,1);
% y1=sol(2,1,1);
% y2=sol(3,1,1);
% y3=sol(4,1,1);
% y4=sol(5,1,1);

c0=0.4;
ini=c0;
for j=1:5000
    A=[1 (-(d0-c0))^(1/2) (-(d0-c0)) (-(d0-c0))^(3/2);...
       1 (-(d1-c0))^(1/2) (-(d1-c0)) (-(d1-c0))^(3/2);...
       1 (-(d2-c0))^(1/2) (-(d2-c0)) (-(d2-c0))^(3/2);...
       1 (-(d3-c0))^(1/2) (-(d3-c0)) (-(d3-c0))^(3/2)];
   m=[y0 y1 y2 y3]';
   coeff=A\m;
   a0=coeff(1);
   a1=coeff(2);
   a2=coeff(3);
   a3=coeff(4);
   %the format a0+a1(d0-c0)^(1/2)... 
   %c0=c0-(a0+a1*(d4-c0)^(1/2)+a2*(d4-c0)+a3*(d4-c0)^(3/2)-y4)/(-a1/2*(d4-c0)^(-1/2)-a2-3/2*a3*(d4-c0)^(1/2));
   
   %the format a0+a1(c0-d0)^(1/2)... 
   c0=c0-(a0+a1*(-(d4-c0))^(1/2)+a2*(-(d4-c0))+a3*(-(d4-c0))^(3/2)-y4)/(a1/2*(-(d4-c0))^(-1/2)+a2+3/2*a3*(c0-d4)^(1/2));
   
   if c0<d4
       j
       break;
   end
end
%%
delta1=-0.0:0.001:0.5;
c=c0;
for i=1:length(delta1)
    %m(i)=a0+a1*(delta1(i)-c)^(1/2)+a2*(delta1(i)-c)+a3*(delta1(i)-c)^(3/2);
    m(i)=a0+a1*(-(delta1(i)-c))^(1/2)+a2*(-(delta1(i)-c))+a3*(-(delta1(i)-c))^(3/2);
end
%=====For Exp1====
% load('Ex1_Sol.mat')
% n=length(delta1);
% n=n+7;
% 
% figure
% plot(delta1,m,'LineWidth',2)
% hold on
% plot(delta(8:n),sold(8:n,1),'LineWidth',2)
% hold on
% plot(delta(8),sold(8,1),'r*','LineWidth',4)

%========For Exp 2=====
figure

plot(delta1,real(m),'LineWidth',2)
hold on
delta_data=-0.05:0.005:0.05;
plot(delta_data,sol(1:end,ind2,1),'LineWidth',2)
hold off
legend('Approx','Exact')
title(['Inital c_0:',num2str(ini)])
figure
plot(delta1,imag(m),'LineWidth',2)
ylabel('Imag(m)')