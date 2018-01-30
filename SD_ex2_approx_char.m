% This code is used to approximate the solution to the 2D SDP along the
% charateristic lines.

clearvars -except x 
format long
load('Exp2_Jacob.mat')

load('Ex2_sol.mat')
% assign the values along line x=y
% delta1=-0.05:0.005:0.05;
% ind1=17;ind2=1;
% dist=delta1(ind2)-delta1(ind1);
% Start=delta1(ind1);


Start=0.0;
xcon=0.0;
%x=Exp2_solve(0.0,-0.12);
%x=sol(21,1,:);
dist=-0.12;
End=0.91;
del_v=0.001;    
del_app=Start:del_v:End;
L2=1;
alp=0;
beta=1;
L1=floor((End-del_app(1))/del_v+1);
m_n=zeros(L1,L2);
y_1=zeros(L1,L2);
y_2=zeros(L1,L2);
x_11=zeros(L1,L2);
x_12=zeros(L1,L2);
x_13=zeros(L1,L2);
x_22=zeros(L1,L2);
x_23=zeros(L1,L2);
x_33=zeros(L1,L2);
s_11=zeros(L1,L2);
s_12=zeros(L1,L2);
s_13=zeros(L1,L2);
s_22=zeros(L1,L2);
s_23=zeros(L1,L2);
s_33=zeros(L1,L2);


% m_n(1) = sol(ind1,ind2,1);
% y_1(1) = sol(ind1,ind2,2);
% y_2(1) = sol(ind1,ind2,3);
% s_11(1) = sol(ind1,ind2,4);
% s_12(1) = sol(ind1,ind2,5);
% s_13(1) = sol(ind1,ind2,6);
% s_22(1) = sol(ind1,ind2,7);
% s_23(1) = sol(ind1,ind2,8);
% s_33(1) = sol(ind1,ind2,9);
% x_11(1) = sol(ind1,ind2,10);
% x_12(1) = sol(ind1,ind2,11);
% x_13(1) = sol(ind1,ind2,12);
% x_22(1) = sol(ind1,ind2,13);
% x_23(1) = sol(ind1,ind2,14);
% x_33(1) = sol(ind1,ind2,15);


m_n(1) = x(1);
y_1(1) = x(2);
y_2(1) = x(3);
s_11(1) = x(4);
s_12(1) = x(5);
s_13(1) = x(6);
s_22(1) = x(7);
s_23(1) = x(8);
s_33(1) = x(9);
x_11(1) = x(10);
x_12(1) = x(11);
x_13(1) = x(12);
x_22(1) = x(13);
x_23(1) = x(14);
x_33(1) = x(15);
for i=1:L1-1
    del_app(i+1)=del_app(1)+del_v*i;
    clear m y1 y2 s1* s2* s3* x1* x2* x3*
    y1=y_1(i);
    y2=y_2(i);
    s11=s_11(i);
    s12=s_12(i);
    s13=s_13(i);
    s22=s_22(i);
    s23=s_23(i);
    s33=s_33(i);
    x11=x_11(i);
    x12=x_12(i);
    x13=x_13(i);
    x22=x_22(i);
    x23=x_23(i);
    x33=x_33(i);
    
    %%=====Test with true derivative===============
    %         y1=sol(i+ind1-1,1+ind2-1,2);
    %         y2=sol(i+ind1-1,1+ind2-1,3);
    %         s11=sol(i+ind1-1,1+ind2-1,4);
    %         s12=sol(i+ind1-1,1+ind2-1,5);
    %         s13=sol(i+ind1-1,1+ind2-1,6);
    %         s22=sol(i+ind1-1,1+ind2-1,7);
    %         s23=sol(i+ind1-1,1+ind2-1,8);
    %         s33=sol(i+ind1-1,1+ind2-1,9);
    %         x11=sol(i+ind1-1,1+ind2-1,10);
    %         x12=sol(i+ind1-1,1+ind2-1,11);
    %         x13=sol(i+ind1-1,1+ind2-1,12);
    %         x22=sol(i+ind1-1,1+ind2-1,13);
    %         x23=sol(i+ind1-1,1+ind2-1,14);
    %         x33=sol(i+ind1-1,1+ind2-1,15);
    %============================================
%     del1=del_app(i);
%     del2=alp*del1+dist;
    del1=beta*del_app(i)+xcon;
    del2=alp*del_app(i)+dist;
    
    J=eval(Jacob_F_V);
    Der_1=eval(Der_F_del1); 
    Der_2=eval(Der_F_del2);

    G=-beta*(J\Der_1')-alp*(J\Der_2');

    
    m_n(i+1)=m_n(i)+del_v*G(1);
    y_1(i+1)=y_1(i)+del_v*G(2);
    y_2(i+1)=y_2(i)+del_v*G(3);
    s_11(i+1)=s_11(i)+del_v*G(4);
    s_12(i+1)=s_12(i)+del_v*G(5);
    s_13(i+1)=s_13(i)+del_v*G(6);
    s_22(i+1)=s_22(i)+del_v*G(7);
    s_23(i+1)=s_23(i)+del_v*G(8);
    s_33(i+1)=s_33(i)+del_v*G(9);
   
    x_11(i+1)=x_11(i)+del_v*G(10);
    x_12(i+1)=x_12(i)+del_v*G(11);
    x_13(i+1)=x_13(i)+del_v*G(12);
    x_22(i+1)=x_22(i)+del_v*G(13);
    x_23(i+1)=x_23(i)+del_v*G(14);
    x_33(i+1)=x_33(i)+del_v*G(15);
    
    clear m y1 y2 s1* s2* s3* x1* x2* x3*
    m=m_n(i+1);
    y1=y_1(i+1);
    y2=y_2(i+1);
    s11=s_11(i+1);
    s12=s_12(i+1);
    s13=s_13(i+1);
    s22=s_22(i+1);
    s23=s_23(i+1);
    s33=s_33(i+1);
    x11=x_11(i+1);
    x12=x_12(i+1);
    x13=x_13(i+1);
    x22=x_22(i+1);
    x23=x_23(i+1);
    x33=x_33(i+1);
    %%=====Test with true derivative===============
    %     m=sol(i+ind1,1+ind2,1);
    %     y1=sol(i+ind1,1+ind2,2);
    %     y2=sol(i+ind1,1+ind2,3);
    %     s11=sol(i+ind1,1+ind2,4);
    %     s12=sol(i+ind1,1+ind2,5);
    %     s13=sol(i+ind1,1+ind2,6);
    %     s22=sol(i+ind1,1+ind2,7);
    %     s23=sol(i+ind1,1+ind2,8);
    %     s33=sol(i+ind1,1+ind2,9);
    %     x11=sol(i+ind1,1+ind2,10);
    %     x12=sol(i+ind1,1+ind2,11);
    %     x13=sol(i+ind1,1+ind2,12);
    %     x22=sol(i+ind1,1+ind2,13);
    %     x23=sol(i+ind1,1+ind2,14);
    %     x33=sol(i+ind1,1+ind2,15);
    %============================================
    del1=beta*del_app(i+1)+xcon;
    del2=alp*del_app(i+1)+dist;
    J=eval(Jacob_F_V) ;
    Der_1=eval(Der_F_del1); Der_2=eval(Der_F_del2);
    G1=-beta*(J\Der_1')-alp*(J\Der_2');
    
    Err_0(i)=eval(f0);
    Err_1(i)=eval(f1);
    Err_2(i)=eval(f2);
    Err_3(i)=eval(f3);
    Err_4(i)=eval(f4);
    Err_5(i)=eval(f5);
    Err_6(i)=eval(f6);
    Err_7(i)=eval(f7);
    Err_8(i)=eval(f8);
    Err_9(i)=eval(f9);
    Err_10(i)=eval(f10);
    Err_11(i)=eval(f11);
    Err_12(i)=eval(f12);
    Err_13(i)=eval(f13);
    Err_14(i)=eval(f14);
    m_n(i+1)=m_n(i)+0.5*del_v*(G(1)+G1(1));
    y_1(i+1)=y_1(i)+0.5*del_v*(G(2)+G1(2));
    y_2(i+1)=y_2(i)+0.5*del_v*(G(3)+G1(3));
    s_11(i+1)=s_11(i)+0.5*del_v*(G(4)+G1(4));
    s_12(i+1)=s_12(i)+0.5*del_v*(G(5)+G1(5));
    s_13(i+1)=s_13(i)+0.5*del_v*(G(6)+G1(6));
    s_22(i+1)=s_22(i)+0.5*del_v*(G(7)+G1(7));
    s_23(i+1)=s_23(i)+0.5*del_v*(G(8)+G1(8));
    s_33(i+1)=s_33(i)+0.5*del_v*(G(9)+G1(9));
    x_11(i+1)=x_11(i)+0.5*del_v*(G(10)+G1(10));
    x_12(i+1)=x_12(i)+0.5*del_v*(G(11)+G1(11));
    x_13(i+1)=x_13(i)+0.5*del_v*(G(12)+G1(12));
    x_22(i+1)=x_22(i)+0.5*del_v*(G(13)+G1(13));
    x_23(i+1)=x_23(i)+0.5*del_v*(G(14)+G1(14));
    x_33(i+1)=x_33(i)+0.5*del_v*(G(15)+G1(15));
    X=real([x_11(i+1) x_12(i+1) x_13(i+1);...
        x_12(i+1) x_22(i+1) x_23(i+1);...
        x_13(i+1) x_23(i+1) x_33(i+1)]);
    S=real([s_11(i+1) s_12(i+1) s_13(i+1);...
        s_12(i+1) s_22(i+1) s_23(i+1);...
        s_13(i+1) s_23(i+1) s_33(i+1)]);
    X_min(i)=min(eig(X));
    S_min(i)=min(eig(S));
         if ((X_min(i))<-10^(-4)) || ((S_min(i))<-10^(-4))
             index = i
             eig_v=[X_min S_min];
             break
         end
    
end


%%
%%=======Define the zone for ploting=======
if (exist('index','var'))
    ind=index;
else
    ind=L1;
end
%%=========================================


figure(15)
plot3(beta*del_app(1:ind)+xcon,alp*del_app(1:ind)+dist,m_n(1:ind),'LineWidth',2)
%plot(beta*del_app(1:ind)+xcon,m_n(1:ind),'LineWidth',2)

hold on
clear m_t
N=floor((0.1-abs(dist))/0.005)+1;
% for i=1:N
%     if ind2>ind1
%         m_t(i)=sol(i+0,i+(ind2-ind1),1);
%        
%     else 
%         m_t(i)=sol(i+(ind1-ind2),i+0,1);
%     end
% end
% inteval=max(0,ind1-ind2);
% delta2=delta1(1+inteval):0.005:delta1(inteval+N);
% plot(delta2,m_t,'--','LineWidt    h',2)
% % hold off
% legend('Approx','Exact')
xlabel(['\delta_1=\delta_2-(',num2str(dist),')'])