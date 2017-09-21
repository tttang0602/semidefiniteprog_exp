syms del s11 s12 s22 x11 x12 x22 y m
load('Ex1_Sol.mat')
sol_d=sold(any(sold,2),:);
ind=any(sold,2);
delta=delta(ind);
delta_v=delta(2)-delta(1);
 n=length(delta);
  dm_ddelta=zeros(n,1);
  %%
  %%====Calculate the derivative using the Davidenko differential equation
 for i =1:n
     del=delta(i);
     m=sol_d(i,1);
     y=sol_d(i,2);
     s11=sol_d(i,3);
     s12=sol_d(i,4);
     s22=sol_d(i,5);
     x11=sol_d(i,6);
     x12=sol_d(i,7);
     x22=sol_d(i,8);
    Jacob_F_V= [1 -1 0 0 0 0 0 0;...
                0  0 0 0 0 -2-del -6-4*del 1-3*del;...
                0 2+del -1 0 0 0 0 0;...
                0 3+2*del 0 -1 0 0 0 0;...
                0 -1+3*del 0 0 -1 0 0 0;...
                0 0 x11 x12 0 s11 s12 0;...
                0 0 x12 x22 0 0 s11 s12;...
                0 0 0 x12 x22 0 s12 s22];

    Jacob_F_del=[0;x11-4*x12-3*x22; y+2; 2*y-1; 3*y+3; 0; 0; 0];        
    %Inv_1=inv(Jacob_F_V);
    %dm=-Inv_1*Jacob_F_del;
    dm=Jacob_F_V\Jacob_F_del;
    dm_ddelta(i)=double(dm(1));
    dy_ddelta(i)=double(dm(2));
    ds11_ddel(i)=double(dm(3));
    ds12_ddel(i)=double(dm(4));
    ds22_ddel(i)=double(dm(5));
    dx11_ddel(i)=double(dm(6));
    dx12_ddel(i)=double(dm(7));
    dx22_ddel(i)=double(dm(8));
 end
 
 
 m_appx(1)=sol_d(end,1);
 y_appx(1)=sol_d(end,2);
 s11_appx(1)=sol_d(end,3);
  s12_appx(1)=sol_d(end,4);
 s22_appx(1)=sol_d(end,5);
 x11_appx(1)=sol_d(end,6);
 x12_appx(1)=sol_d(end,7);
 x22_appx(1)=sol_d(end,8);
%% 
save('Exl_Der.mat','*del*')
 
 %%
clear all 

 syms del s11 s12 s22 x11 x12 x22 y m
load('Ex1_Sol.mat')
load('Exl_Der.mat')
sol_d=sold(any(sold,2),:);
%ind=any(sold,2);
%delta1=delta(ind);
startn=1000;
  m_appx1(1)=sol_d(startn,1);
  y_appx1(1)=sol_d(startn,2);
 s11_appx1(1)=sol_d(startn,3);
  s12_appx1(1)=sol_d(startn,4);
 s22_appx1(1)=sol_d(startn,5);
 x11_appx1(1)=sol_d(startn,6);
 x12_appx1(1)=sol_d(startn,7);
 x22_appx1(1)=sol_d(startn,8);
 dm_ddelta1(1)=dm_ddelta(startn);
 dy_ddelta1(i)=dy_ddelta(startn);
 
 ds111_ddel(1)=ds11_ddel(startn);
 ds121_ddel(1)=ds12_ddel(startn);
 ds221_ddel(1)=ds22_ddel(startn);
 dx111_ddel(1)=dx11_ddel(startn);
 dx121_ddel(1)=dx12_ddel(startn);
 dx221_ddel(1)=dx22_ddel(startn);

 k=startn*10+150;
 %del_app=zeros(k,1);
 del_app(1)=delta(startn);
 del_v=0.0001;
 %%
 for i=1:k
     del_app(i+1)=delta(startn)-i*del_v;
     m_appx1(i+1)=m_appx1(i)-del_v*dm_ddelta1(i);
     y_appx1(i+1)=y_appx1(i)-del_v*dy_ddelta1(i);
     s11_appx1(i+1)=s11_appx1(i)-del_v*ds111_ddel(i);
     s12_appx1(i+1)=s12_appx1(i)-del_v*ds121_ddel(i);
     s22_appx1(i+1)=s22_appx1(i)-del_v*ds221_ddel(i);
     x11_appx1(i+1)=x11_appx1(i)-del_v*dx111_ddel(i);
     x12_appx1(i+1)=x12_appx1(i)-del_v*dx121_ddel(i);
     x22_appx1(i+1)=x22_appx1(i)-del_v*dx221_ddel(i);
     
     
     
     del=del_app(i+1);
     m=m_appx1(i+1);
     y=y_appx1(i+1);
     s11=s11_appx1(i+1);
     s12=s12_appx1(i+1);
     s22=s22_appx1(i+1);
     x11=x11_appx1(i+1);
     x12=x12_appx1(i+1);
     x22=x22_appx1(i+1);
    Jacob_F_V= [1 -1 0 0 0 0 0 0;...
                0  0 0 0 0 -2-del -6-4*del 1-3*del;...
                0 2+del -1 0 0 0 0 0;...
                0 3+2*del 0 -1 0 0 0 0;...
                0 -1+3*del 0 0 -1 0 0 0;...
                0 0 x11 x12 0 s11 s12 0;...
                0 0 x12 x22 0 0 s11 s12;...
                0 0 0 x12 x22 0 s12 s22];

    Jacob_F_del=[0;x11-4*x12-3*x22; y+2; 2*y-1; 3*y+3; 0; 0; 0];        
    dm=-Jacob_F_V\Jacob_F_del;
    dm_ddelta1(i+1)=double(dm(1));
    dy_ddelta1(i+1)=double(dm(2));
    ds111_ddel(i+1)=double(dm(3));
    ds121_ddel(i+1)=double(dm(4));
    ds221_ddel(i+1)=double(dm(5));
    dx111_ddel(i+1)=double(dm(6));
    dx121_ddel(i+1)=double(dm(7));
    dx221_ddel(i+1)=double(dm(8));
    
    m_appx1(i+1)=m_appx1(i)-0.5*del_v*(dm_ddelta1(i)+dm_ddelta1(i+1));
    y_appx1(i+1)=y_appx1(i)-0.5*del_v*(dy_ddelta1(i)+dy_ddelta1(i+1));
    
    s11_appx1(i+1)=s11_appx1(i)-0.5*del_v*(ds111_ddel(i)+ds111_ddel(i+1));
     s12_appx1(i+1)=s12_appx1(i)-0.5*del_v*(ds121_ddel(i)+ds121_ddel(i+1));
     s22_appx1(i+1)=s22_appx1(i)-0.5*del_v*(ds221_ddel(i)+ds221_ddel(i+1));
     x11_appx1(i+1)=x11_appx1(i)-0.5*del_v*(dx111_ddel(i)+dx111_ddel(i+1));
     x12_appx1(i+1)=x12_appx1(i)-0.5*del_v*(dx121_ddel(i)+dx121_ddel(i+1));
     x22_appx1(i+1)=x22_appx1(i)-0.5*del_v*(dx221_ddel(i)+dx221_ddel(i+1));
     
     S=[s11_appx1(i+1) s12_appx1(i+1);s12_appx1(i+1) s22_appx1(i+1)];
     X=[x11_appx1(i+1) x12_appx1(i+1);x12_appx1(i+1) x22_appx1(i+1)];
     S_min=min(eig(S));
     X_min=min(eig(X));
     Err(i)=x11_appx1(i+1)*s11_appx1(i+1)+x12_appx1(i+1)*s12_appx1(i+1);
     Err1(i)=1+2*y_appx1(i+1)-s11_appx1(i+1)+del*(2+y_appx1(i+1));
     
     if ((X_min)<-10^(-5)) || ((S_min)<-10^(-5))
         i
         break
     end
 end
 %%

 figure 
 plot(del_app,y_appx1,'LineWidth',2)
 hold on
 plot(delta(1:startn),sol_d(1:startn,2),'LineWidth',2)
 legend('Appriximation','Feasible solution')
 title(['\delta_0=',num2str(delta(startn))])  
 
 
 
 %%
%%===Approximation using Puiseux series===

%%==Using 4 points mesh size 0.001
clear all

 d0=-0.489;
d1=-0.49;
d2=-0.491;
d3=-0.492;
%d4=-0.493;
y0=0.204305518472662;
y1=0.201107684592419;
y2=0.197482095669551;
y3=0.193115453265627;
%y4=0.186900321003191;
%%Using Only 4 points..===

a0=1.798669263570130e-01 ;
a1=3.320255563575221e-01;
a2=5.388985688134520e-01 ;
c=-4.934142632055215e-01;
delta1=-0.493:0.001:0;
n=length(delta1);
for i=1:length(delta1)
    m(i)=a0+a1*(delta1(i)-c)^(1/2)+a2*(delta1(i)-c);
end
load('Ex1_Sol.mat')
figure
plot(delta1,m,'LineWidth',2)
hold on

n=n+7;
plot(delta(8:n),sold(8:n,1),'LineWidth',2)
hold on
plot(delta(8),sold(8,1),'r*','LineWidth',4)

legend('Puiseux Appriximation','Feasible solution')

%%
%%Using mesh size 0.01 and 5 points
a0=1.802679588485794e-01 ;
a1=3.375148085762173e-01 ;
a2=3.676669254170282e-01 ;
a3=1.050133576182947e+00 ;
c=-4.933113752888860e-01 ;

for i=1:length(delta1)
    m(i)=a0+a1*(delta1(i)-c)^(1/2)+a2*(delta1(i)-c)*2+a3*(delta1(i)-c)^(3/2)*6;
end
load('Ex1_Sol.mat')
ind=find(delta==0);
figure
plot(delta1,m,'LineWidth',2)
hold on
n=length(delta1);

n=n+7;
plot(delta(8:n),sold(8:n,1),'LineWidth',2)
hold on
plot(delta(8),sold(8,1),'r*','LineWidth',4)

legend('Puiseux Appriximation','Feasible solution')

%%
%%==Using mesh size 0.1 and 5 points
clear m
a0=2.052792779786183e-01;
a1=1.329096303493892e-01 ;
a2=9.841560034079134e-01 ;
a3=3.899223787794775e-01 ;
c=-4.923425657594596e-01;
t0=-4.840005949851205e-02 ;
t1=3.199102464121139e-01 ;
t2=4.498250390534742e-01 ;
t3=5.498568593365545e-01 ;
t4=6.343047893240754e-01 ;
for i=1:length(delta1)
    m(i)=a0+a1*(delta1(i)-c)^(1/2)+a2*(delta1(i)-c)+a3*(delta1(i)-c)^(3/2);
end
load('Ex1_Sol.mat')
figure
plot(delta1,m,'LineWidth',2)
hold on
n=length(delta1);
n=n+7;
plot(delta(8:n),sold(8:n,1),'LineWidth',2)
hold on
plot(delta(8),sold(8,1),'r*','LineWidth',4)

legend('Puiseux Appriximation','Feasible solution')
%%
%%==Using mesh size 0.01 and 5 points (away from singularity)
clear m
a0=0.196047318492610154827e0;
a1=0.248677121455738156655e0;
a2=0.678086673104939838307e0;
a3=0.661252751531543495956e0;
c=-0.488470114307021013064e0;
t0=0.297439261542623132391e0;
t1=0.280125176139205593812e0;
t2=0.261667946655720445284e0;
t3=0.241805943489859155878e0;
t4=0.220159293028981548274e0;
t0^2+(0.4)+c
t1^2+0.41+c
t2^2+0.42+c
t3^2+0.43+c
for i=1:length(delta1)
    m(i)=a0+a1*(delta1(i)-c)^(1/2)+a2*(delta1(i)-c)+a3*(delta1(i)-c)^(3/2);
end
load('Ex1_Sol.mat')
n=length(delta1);
n=n+7;
figure
plot(delta1(6:end),m(6:end),'LineWidth',2)
hold on
plot(delta(8:n),sold(8:n,1),'LineWidth',2)
hold on
plot(delta1(6),m(6),'m*','LineWidth',4)
hold on
plot(delta(8),sold(8,1),'r*','LineWidth',4)

legend('Puiseux Appriximation','Feasible solution')
%%
%%===Using 5 points mesh size 0.001====
a0 = 0.179076316945588857871e0 ;
a1 = 0.354246832213496549047e0;
a2 = 0.255839250279742218497;
a3 = 0.133054991274314350679e1 ;
c = -0.493471218563230327602e0;

delta1=-0.493:0.001:0;

for i=1:length(delta1)
    m(i)=a0+a1*(delta1(i)-c)^(1/2)+a2*(delta1(i)-c)+a3*(delta1(i)-c)^(3/2);
end
load('Ex1_Sol.mat')
n=length(delta1);
n=n+7;

figure
plot(delta1,m,'LineWidth',2)
hold on
plot(delta(8:n),sold(8:n,1),'LineWidth',2)
hold on
plot(delta(8),sold(8,1),'r*','LineWidth',4)

legend('Puiseux Appriximation','Feasible solution')

%%
a0 = 1.960473184919618e-01;
a1 = 2.486771214580826e-01;
a2 = 6.780866730986450e-01 ;    
a3 = 6.612527515372117e-01 ;
c = -4.884701143072788e-01 ;
t0 = 2.974392615430565e-01 ;
t1 = 2.801251761396658e-01 ;
t2 = 2.616679466562131e-01 :
t3 = 2.418059434903923e-01 :
t4 = 2.201592930295671e-01 ;
