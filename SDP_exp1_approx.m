%===========This code using modified eulur method to approximate 1D
%solution and plot results of Puiseux seris apprixmation.


%%%%========================================================================
%%%%%================Apprixmation
%%%%%section==============================================================

clear all
syms del s11 s12 s22 x11 x12 x22 y m

load('Ex1_Sol.mat')
load('Der_Ex1.mat')
sol_d=sold(any(sold,2),:);
mexact=@(d) (13*d^2+sqrt(189*d^4+654*d^3+959*d^2+590*d+125)+19*d+9)/(2*(d^2+7*d+11));
%==============Positive+==========
% del_v=0.005;
% End=2.25;
% del_app=0:del_v:End;
% L=length(del_app);
%================================

%===========Negative================
del_v=-0.0001;
End=-0.6;
del_app=0:del_v:End;
L=length(del_app);
%===================================

%Initialize solution value
m_appx1=zeros(L,1);
y_appx1=zeros(L,1);
s11_appx1=zeros(L,1);
s12_appx1=zeros(L,1);
s22_appx1=zeros(L,1);
x11_appx1=zeros(L,1);
x12_appx1=zeros(L,1);
x22_appx1=zeros(L,1);

%Initialize derivatives value
dm_ddelta1=zeros(L,1);
dy_ddelta1=zeros(L,1);
ds111_ddel=zeros(L,1);
ds121_ddel=zeros(L,1);
ds221_ddel=zeros(L,1);
dx111_ddel=zeros(L,1);
dx121_ddel=zeros(L,1);
dx221_ddel=zeros(L,1);

%Assign boundary value to solution
startn=501;
  m_appx1(1)=sold(startn,1);
  y_appx1(1)=sold(startn,2);
 s11_appx1(1)=sold(startn,3);
  s12_appx1(1)=sold(startn,4);
 s22_appx1(1)=sold(startn,5);
 x11_appx1(1)=sold(startn,6);
 x12_appx1(1)=sold(startn,7);
 x22_appx1(1)=sold(startn,8);
 
 %Assign value to derivative of the solution boundary value
 m = m_appx1(1);
 y = y_appx1(1);
 s11 = s11_appx1(1);
 s12 = s12_appx1(1);
 s22 = s22_appx1(1);
 x11 = x11_appx1(1);
 x12 = x12_appx1(1);
 x22 = x22_appx1(1);
 Jacob_F_V= [1 -1 0 0 0 0 0 0;...
                0  0 0 0 0 -2 -6 1;...
                0 2 -1 0 0 0 0 0;...
                0 3 0 -1 0 0 0 0;...
                0 -1 0 0 -1 0 0 0;...
                0 0 x11 x12 0 s11 s12 0;...
                0 0 x12 x22 0 0 s11 s12;...
                0 0 0 x12 x22 0 s12 s22];

Jacob_F_del=[0;x11-4*x12-3*x22; y+2; 2*y-1; 3*y+3; 0; 0; 0];
Dm=-Jacob_F_V\Jacob_F_del;

dm_ddelta1(1)=Dm(1);
dy_ddelta1(1)=Dm(2);
 ds111_ddel(1)=Dm(3);
 ds121_ddel(1)=Dm(4);
 ds221_ddel(1)=Dm(5);
 dx111_ddel(1)=Dm(6);
 dx121_ddel(1)=Dm(7);
 dx221_ddel(1)=Dm(8);
    dlambdax=zeros(L,1);
    dlambdas=zeros(L,1);
 %%
 for i=1:L-1
     
     m_exp(i)=(del_app(i)^2*13+19*del_app(i)+9+sqrt(189*del_app(i)^4+654*del_app(i)^3+959*del_app(i)^2+590*del_app(i)+125))/(2*del_app(i)^2+14*del_app(i)+22);
     
     m_appx1(i+1)=m_appx1(i)+del_v*dm_ddelta1(i);
     y_appx1(i+1)=y_appx1(i)+del_v*dy_ddelta1(i);
     s11_appx1(i+1)=s11_appx1(i)+del_v*ds111_ddel(i);
     s12_appx1(i+1)=s12_appx1(i)+del_v*ds121_ddel(i);
     s22_appx1(i+1)=s22_appx1(i)+del_v*ds221_ddel(i);
     x11_appx1(i+1)=x11_appx1(i)+del_v*dx111_ddel(i);
     x12_appx1(i+1)=x12_appx1(i)+del_v*dx121_ddel(i);
     x22_appx1(i+1)=x22_appx1(i)+del_v*dx221_ddel(i);
     
     
     
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
    Dm=-Jacob_F_V\Jacob_F_del;
    dm_ddelta1(i+1)=double(Dm(1));
    dy_ddelta1(i+1)=double(Dm(2));
    ds111_ddel(i+1)=double(Dm(3));
    ds121_ddel(i+1)=double(Dm(4));
    ds221_ddel(i+1)=double(Dm(5));
    dx111_ddel(i+1)=double(Dm(6));
    dx121_ddel(i+1)=double(Dm(7));
    dx221_ddel(i+1)=double(Dm(8));
    
    m_appx1(i+1)=m_appx1(i)+0.5*del_v*(dm_ddelta1(i)+dm_ddelta1(i+1));
    y_appx1(i+1)=y_appx1(i)+0.5*del_v*(dy_ddelta1(i)+dy_ddelta1(i+1));
    
    s11_appx1(i+1)=s11_appx1(i)+0.5*del_v*(ds111_ddel(i)+ds111_ddel(i+1));
     s12_appx1(i+1)=s12_appx1(i)+0.5*del_v*(ds121_ddel(i)+ds121_ddel(i+1));
     s22_appx1(i+1)=s22_appx1(i)+0.5*del_v*(ds221_ddel(i)+ds221_ddel(i+1));
     x11_appx1(i+1)=x11_appx1(i)+0.5*del_v*(dx111_ddel(i)+dx111_ddel(i+1));
     x12_appx1(i+1)=x12_appx1(i)+0.5*del_v*(dx121_ddel(i)+dx121_ddel(i+1));
     x22_appx1(i+1)=x22_appx1(i)+0.5*del_v*(dx221_ddel(i)+dx221_ddel(i+1));
%==========================Criterial One================SDP condition===========     
     S=[s11_appx1(i+1) s12_appx1(i+1);s12_appx1(i+1) s22_appx1(i+1)];
     X=[x11_appx1(i+1) x12_appx1(i+1);x12_appx1(i+1) x22_appx1(i+1)];
     Err(i)=x11_appx1(i+1)*s11_appx1(i+1)+x12_appx1(i+1)*s12_appx1(i+1);
     Err1(i)=1+2*y_appx1(i+1)-s11_appx1(i+1)+del*(2+y_appx1(i+1));
     

%======================================================================

        [Vx,Dx,Wx]=eig(X);
        [Vs,Ds,Ws]=eig(S);
        [X_min(i),ox]=min(diag(Dx));
        [S_min(i),os]=min(diag(Ds));
        
%         %==================Criteria eigenvalue=========================
%         if ((X_min(i))<-10^(1))|| ((S_min(i))<-10^(1))
%             i
%             break
%         end
        vx=Vx(:,ox);
        vs=Vs(:,os);
        wx=Wx(:,ox);
        ws=Ws(:,os);
        m=m_appx1(i+1);
        y=y_appx1(i+1);
        s11=s11_appx1(i+1);
        s12=s12_appx1(i+1);
        s22=s22_appx1(i+1);
        x11=x11_appx1(i+1);
        x12=x12_appx1(i+1);
        x22=x22_appx1(i+1);
        G=eval(dm);
        dx=[G(6),G(7);...
            G(7),G(8)];
        ds=[G(3),G(4);...
            G(4),G(5)];
        for p=1:2
            for q=1:2
                dlambdax(i)=dlambdax(i)+wx(p)*vx(q)/(vx'*wx)*dx(p,q);
                dlambdas(i)=dlambdas(i)+ws(p)*vs(q)/(vs'*ws)*ds(p,q);
            end
        end
%=====================Criteria Zero======Derivative of lambda==========        
        if abs(dlambdas(i))>1||abs(dlambdax(i))>1
            index=i
            break
        end
%=============Criterial Two======================================
%==Ratio of second derivative to first derivative================
%  Using exact derivative ======================================
%     del=del_app(i+1);
%      m=m_appx1(i+1);
%      y=y_appx1(i+1);
%      s11=s11_appx1(i+1);
%      s12=s12_appx1(i+1);
%      s22=s22_appx1(i+1);
%      x11=x11_appx1(i+1);
%      x12=x12_appx1(i+1);
%      x22=x22_appx1(i+1);
%      Der1_m=eval(dm(1));
%      Der2_m=eval(dm1(1));

%===Using approximated derivative ============================
%     if i>=2
%         Der1_m=(m_appx1(i+1)-m_appx1(i))/del_v;
%         Der2_m=(m_appx1(i+1)-2*m_appx1(i)+m_appx1(i-1))/(del_v)^2;
%         rat=Der2_m/Der1_m;
%         if abs(rat)>=50
%             i
%             break
%         end
%     end
 end
 %%

 figure 
 plot(delta(8:501),sold(8:501,1),'g-','LineWidth',2)
 hold on
 plot(del_app(1:index),m_exp(1:index),'r-.','LineWidth',2)
 hold on
 plot(del_app(1:index+1),m_appx1(1:index+1),'b--','LineWidth',2)
 hold on
% scatter(-0.49347124103032225,mm)
 legend('Method 1','Method 2','Method 3','Singularity')
 xlim([-0.6 0])
 %title(['\delta_0=',num2str(delta(startn))])  
 xlabel('\delta')
 ylabel('m(\delta)')
 save('Ex1_sol_app_p','del_app','*_appx1','i')
 
 
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


