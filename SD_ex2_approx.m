A1 = [2 4 -9; 4 10 -6; -9 -6 8];
A2 = [4 -6 4; -6 0 -7; 4 -7 4];
C  = [-8 2 8; 2 2 6; 8 6 0];
D1 = [4 5 2; 5 -8 0; 2 0 -6];
E11= [8 6 4 ; 6 -4 2; 4 2 2];
E12= [4 -2 -5; -2 -6 9; -5 9 0];
D2 = [6 -5 6; -5 -6 -3; 6 -3 -4];
E21= [-10 -2 1; -2 10 -4; 1 -4 -10];
E22= [10 -9 2; -9 10 -2; 2 -2 8];

syms m y1 y2 s11 s12 s13 s22 s23 s33 x11 x12 x13 x22 x23 x33
X_sym = [x11 x12 x13;...
         x12 x22 x23;...
         x13 x23 x33];
S_sym = [s11 s12 s13;...
         s12 s22 s23;...
         s13 s23 s33];
ndel=20;
delta1=linspace(-0.05,0.05,ndel);
max_sol=zeros(ndel,ndel);
syms del1 del2
%%===
f0=m-y1-4*y2;
 f1=trace((A1+del1*E11+del2*E21)*X_sym)-1;
 f2=trace((A2+del1*E12+del2*E22)*X_sym)+4;
 matr_stat = C-(y1*A1+y2*A2)+del1*(D1-y1*E11-y2*E12)+del2*(D2-y1*E21-y2*E22)-S_sym;
 f3=matr_stat(1,1);
 f4 = matr_stat(1,2);
f5 = matr_stat(1,3);
f6 = matr_stat(2,2);
f7 = matr_stat(2,3);
f8 = matr_stat(3,3);
eqn01=0;
eqn02=0;
eqn11=diff(f1,del1);
eqn12=diff(f1,del2);
eqn21=diff(f2,del1);
eqn22=diff(f2,del2);
eqn31=diff(f3,del1);
eqn32=diff(f3,del2);
eqn41=diff(f4,del1);
eqn42=diff(f4,del2);

eqn51=diff(f5,del1);
eqn52=diff(f5,del2);

eqn61=diff(f6,del1);
eqn62=diff(f6,del2);

eqn71=diff(f7,del1);
eqn72=diff(f7,del2);

eqn81=diff(f8,del1);
eqn82=diff(f8,del2);

matr_comp=S_sym*X_sym;

f9=matr_comp(1,1);
f10=matr_comp(1,2);
f11=matr_comp(1,3);
f12=matr_comp(2,2);
f13=matr_comp(2,3);
f14=matr_comp(3,3);
eqn91=diff(f9,del1);
eqn92=diff(f9,del2);

eqn101=diff(f10,del1);
eqn102=diff(f10,del2);


eqn111=diff(f11,del1);
eqn112=diff(f11,del2);

eqn121=diff(f12,del1);
eqn122=diff(f12,del2);

eqn131=diff(f13,del1);
eqn132=diff(f13,del2);

eqn141=diff(f14,del1);
eqn142=diff(f14,del2);

%%
Der_F_del1=[eqn01,eqn11,...
    eqn21,eqn31,eqn41,eqn51,eqn61,eqn71,eqn81,...
    eqn91,eqn101,eqn111,eqn121,eqn131,eqn141]';
Der_F_del2=[eqn02,eqn12,...
    eqn22,eqn32,eqn42,eqn52,eqn62,eqn72,eqn82,...
    eqn92,eqn102,eqn112,eqn122,eqn132,eqn142]';

Jacob_F_V=[diff(f0,m),diff(f0,y1),diff(f0,y2),...
            diff(f0,x11),diff(f0,x12),diff(f0,x13),diff(f0,x22),diff(f0,x23),diff(f0,x33),...
            diff(f0,s11),diff(f0,s12),diff(f0,s13),diff(f0,s22),diff(f0,s23),diff(f0,s33);...
           diff(f1,m),diff(f1,y1),diff(f1,y2),...
            diff(f1,x11),diff(f1,x12),diff(f1,x13),diff(f1,x22),diff(f1,x23),diff(f1,x33),...
            diff(f1,s11),diff(f1,s12),diff(f1,s13),diff(f1,s22),diff(f1,s23),diff(f1,s33);...
            diff(f2,m),diff(f2,y1),diff(f2,y2),...
            diff(f2,x11),diff(f2,x12),diff(f2,x13),diff(f2,x22),diff(f2,x23),diff(f2,x33),...
            diff(f2,s11),diff(f2,s12),diff(f2,s13),diff(f2,s22),diff(f2,s23),diff(f2,s33);...
            diff(f3,m),diff(f3,y1),diff(f3,y2),...
            diff(f3,x11),diff(f3,x12),diff(f3,x13),diff(f3,x22),diff(f3,x23),diff(f3,x33),...
            diff(f3,s11),diff(f3,s12),diff(f3,s13),diff(f3,s22),diff(f3,s23),diff(f3,s33);...
            diff(f4,m),diff(f4,y1),diff(f4,y2),...
            diff(f4,x11),diff(f4,x12),diff(f4,x13),diff(f4,x22),diff(f4,x23),diff(f4,x33),...
            diff(f4,s11),diff(f4,s12),diff(f4,s13),diff(f4,s22),diff(f4,s23),diff(f4,s33);...
            diff(f5,m),diff(f5,y1),diff(f5,y2),...
            diff(f5,x11),diff(f5,x12),diff(f5,x13),diff(f5,x22),diff(f5,x23),diff(f5,x33),...
            diff(f5,s11),diff(f5,s12),diff(f5,s13),diff(f5,s22),diff(f5,s23),diff(f5,s33);...
            diff(f6,m),diff(f6,y1),diff(f6,y2),...
            diff(f6,x11),diff(f6,x12),diff(f6,x13),diff(f6,x22),diff(f6,x23),diff(f6,x33),...
            diff(f6,s11),diff(f6,s12),diff(f6,s13),diff(f6,s22),diff(f6,s23),diff(f6,s33);...
            diff(f7,m),diff(f7,y1),diff(f7,y2),...
            diff(f7,x11),diff(f7,x12),diff(f7,x13),diff(f7,x22),diff(f7,x23),diff(f7,x33),...
            diff(f7,s11),diff(f7,s12),diff(f7,s13),diff(f7,s22),diff(f7,s23),diff(f7,s33);...
            diff(f8,m),diff(f8,y1),diff(f8,y2),...
            diff(f8,x11),diff(f8,x12),diff(f8,x13),diff(f8,x22),diff(f8,x23),diff(f8,x33),...
            diff(f8,s11),diff(f8,s12),diff(f8,s13),diff(f8,s22),diff(f8,s23),diff(f8,s33);...
            diff(f9,m),diff(f9,y1),diff(f9,y2),...
            diff(f9,x11),diff(f9,x12),diff(f9,x13),diff(f9,x22),diff(f9,x23),diff(f9,x33),...
            diff(f9,s11),diff(f9,s12),diff(f9,s13),diff(f9,s22),diff(f9,s23),diff(f9,s33);...
            diff(f10,m),diff(f10,y1),diff(f10,y2),...
            diff(f10,x11),diff(f10,x12),diff(f10,x13),diff(f10,x22),diff(f10,x23),diff(f10,x33),...
            diff(f10,s11),diff(f10,s12),diff(f10,s13),diff(f10,s22),diff(f10,s23),diff(f10,s33);...
            diff(f11,m),diff(f11,y1),diff(f11,y2),...
            diff(f11,x11),diff(f11,x12),diff(f11,x13),diff(f11,x22),diff(f11,x23),diff(f11,x33),...
            diff(f11,s11),diff(f11,s12),diff(f11,s13),diff(f11,s22),diff(f11,s23),diff(f11,s33);...
            diff(f12,m),diff(f12,y1),diff(f12,y2),...
            diff(f12,x11),diff(f12,x12),diff(f12,x13),diff(f12,x22),diff(f12,x23),diff(f12,x33),...
            diff(f12,s11),diff(f12,s12),diff(f12,s13),diff(f12,s22),diff(f12,s23),diff(f12,s33);...
            diff(f13,m),diff(f13,y1),diff(f13,y2),...
            diff(f13,x11),diff(f13,x12),diff(f13,x13),diff(f13,x22),diff(f13,x23),diff(f13,x33),...
            diff(f13,s11),diff(f13,s12),diff(f13,s13),diff(f13,s22),diff(f13,s23),diff(f13,s33);...
            diff(f14,m),diff(f14,y1),diff(f14,y2),...
            diff(f14,x11),diff(f14,x12),diff(f14,x13),diff(f14,x22),diff(f14,x23),diff(f14,x33),...
            diff(f14,s11),diff(f14,s12),diff(f14,s13),diff(f14,s22),diff(f14,s23),diff(f14,s33)];

%%
save('Exp2_Jacob.mat','Jacob_F_V','Der_F_del1','Der_F_del2');
%%
clear all
load('Exp2_Jacob.mat')
load('Ex2_sol.mat')
%load('Ex2_sol_fine.mat')
%%

delta1=(-0.05:0.005:0.05);
delta2=(-0.05:0.005:0.05);
del_1=delta1(2)-delta1(1);
del_2=delta2(2)-delta2(1);
L1=length(delta1);
L2=length(delta2);
m_n=zeros(L1,L2);
y_1=zeros(L1,L2);
y_2=zeros(L1,L2);
x_11=zeros(L1,L2);
x_12=zeros(L1,L2);
x_13=zeros(L1,L2);
x_22=zeros(L1,L2);z`
x_23=zeros(L1,L2);
x_33=zeros(L1,L2);
s_11=zeros(L1,L2);
s_12=zeros(L1,L2);
s_13=zeros(L1,L2);
s_22=zeros(L1,L2);
s_23=zeros(L1,L2);
s_33=zeros(L1,L2);
ind1=1;ind2=1;
%Scheme boundary delta1(ind1) with delta2(ind2:end) and delta2(ind1) with
%delta1(ind2:end)
%sol=zeros(L1,L2,15);
%sol(1,:,:)=sol_del2(1,:,:);
%sol(:,1,:)=sol_del1(:,1,:);
m_n(1,:)=sol(ind1,ind2:end,1);
m_n(:,1)=sol(ind2:end,ind1,1);
y_1(1,:)=sol(ind1,ind2:end,2);
y_1(:,1)=sol(ind2:end,ind1,2);
y_2(1,:)=sol(ind1,ind2:end,3);
y_2(:,1)=sol(ind2:end,ind1,3);
s_11(1,:)=sol(ind1,ind2:end,4);
s_11(:,1)=sol(ind2:end,ind1,4);
s_12(1,:)=sol(ind1,ind2:end,5);
s_12(:,1)=sol(ind2:end,ind1,5);
s_13(1,:)=sol(ind1,ind2:end,6);
s_13(:,1)=sol(ind2:end,ind1,6);
s_22(1,:)=sol(ind1,ind2:end,7);
s_22(:,1)=sol(ind2:end,ind1,7);
s_23(1,:)=sol(ind1,ind2:end,8);
s_23(:,1)=sol(ind2:end,ind1,8);
s_33(1,:)=sol(ind1,ind2:end,9);
s_33(:,1)=sol(ind2:end,ind1,9);
x_11(1,:)=sol(ind1,ind2:end,10);
x_11(:,1)=sol(ind2:end,ind1,10);
x_12(1,:)=sol(ind1,ind2:end,11);
x_12(:,1)=sol(ind2:end,ind1,11);
x_13(1,:)=sol(ind1,ind2:end,12);
x_13(:,1)=sol(ind2:end,ind1,12);
x_22(1,:)=sol(ind1,ind2:end,13);
x_22(:,1)=sol(ind2:end,ind1,13);
x_23(1,:)=sol(ind1,ind2:end,14);
x_23(:,1)=sol(ind2:end,ind1,14);
x_33(1,:)=sol(ind1,ind2:end,15);
x_33(:,1)=sol(ind2:end,ind1,15);

% m_n(1,:)=fliplr(m_n(1,:));
% m_n(:,1)=fliplr(m_n(:,1));
% y_1(1,:)=fliplr(y_1(1,:));
% y_1(:,1)=fliplr(y_1(:,1));
% y_2(1,:)=fliplr(y_2(1,:));
% y_2(:,1)=fliplr(y_2(:,1));
% s_11(1,:)=fliplr(s_11(1,:));
% s_11(:,1)=fliplr(s_11(:,1));
% s_12(1,:)=fliplr(s_12(1,:));
% s_12(:,1)=fliplr(s_12(:,1));
% s_13(1,:)=fliplr(s_13(1,:));
% s_13(:,1)=fliplr(s_13(:,1));
% s_22(1,:)=fliplr(s_22(1,:));
% s_22(:,1)=fliplr(s_22(:,1));
% s_23(1,:)=fliplr(s_23(1,:));
% s_23(:,1)=fliplr(s_23(:,1));
% s_33(1,:)=fliplr(s_33(1,:));
% s_33(:,1)=fliplr(s_33(:,1));
% x_11(1,:)=fliplr(x_11(1,:));
% x_11(:,1)=fliplr(x_11(:,1));
% x_12(1,:)=fliplr(x_12(1,:));
% x_12(:,1)=fliplr(x_12(:,1));
% x_13(1,:)=fliplr(x_13(1,:));
% x_13(:,1)=fliplr(x_13(:,1));
% x_22(1,:)=fliplr(x_22(1,:));
% x_22(:,1)=fliplr(x_22(:,1));
% x_23(1,:)=fliplr(x_23(1,:));
% x_23(:,1)=fliplr(x_23(:,1));
% x_33(1,:)=fliplr(x_33(1,:));
% x_33(:,1)=fliplr(x_33(:,1));

hm=zeros(L1,L2);
hy1=zeros(L1,L2);
hy2=zeros(L1,L2);
hx11=zeros(L1,L2);
hx12=zeros(L1,L2);
hx13=zeros(L1,L2);
hx22=zeros(L1,L2);
hx23=zeros(L1,L2);
hx33=zeros(L1,L2);
hs11=zeros(L1,L2);
hs12=zeros(L1,L2);
hs13=zeros(L1,L2);
hs22=zeros(L1,L2);
hs23=zeros(L1,L2);
hs33=zeros(L1,L2);
hm(:,1)=m_n(:,1);
hy1(:,1)=y_1(:,1);
hy2(:,1)=y_2(:,1);
hx11(:,1)=x_11(:,1);
hx12(:,1)=x_12(:,1);
hx13(:,1)=x_13(:,1);
hx22(:,1)=x_22(:,1);
hx23(:,1)=x_23(:,1);
hx33(:,1)=x_33(:,1);
hs11(:,1)=s_11(:,1);
hs12(:,1)=s_12(:,1);
hs13(:,1)=s_13(:,1);
hs22(:,1)=s_22(:,1);
hs23(:,1)=s_23(:,1);
hs33(:,1)=s_33(:,1);

for i=2:L1
    for j=2:L2
        if j<=2
            hm(i-1,j)=m_n(i-1,j);
            hy1(i-1,j)=y_1(i-1,j);
            hy2(i-1,j)=y_2(i-1,j);
            hx11(i-1,j)=x_11(i-1,j);
            hx12(i-1,j)=x_12(i-1,j);
            hx13(i-1,j)=x_13(i-1,j);
            hx22(i-1,j)=x_22(i-1,j);
            hx23(i-1,j)=x_23(i-1,j);
            hx33(i-1,j)=x_33(i-1,j);
        elseif j>=L2-1
            hm(i-1,j)=m_n(i-1,j);
            hy1(i-1,j)=y_1(i-1,j);
            hy2(i-1,j)=y_2(i-1,j);
            hx11(i-1,j)=x_11(i-1,j);
            hx12(i-1,j)=x_12(i-1,j);
            hx13(i-1,j)=x_13(i-1,j);
            hx22(i-1,j)=x_22(i-1,j);
            hx23(i-1,j)=x_23(i-1,j);
            hx33(i-1,j)=x_33(i-1,j);
        else
            hm(i-1,j)=m_n(i-1,j)+1/2*minmod(m_n(i-1,j+1)-m_n(i-1,j),m_n(i-1,j)-m_n(i-1,j-1));
            hy1(i-1,j)=y_1(i-1,j)+1/2*minmod(y_1(i-1,j+1)-y_1(i-1,j),y_1(i-1,j)-y_1(i-1,j-1));
            hy2(i-1,j)=y_2(i-1,j)+1/2*minmod(y_2(i-1,j+1)-y_2(i-1,j),y_2(i-1,j)-y_2(i-1,j-1));
            hx11(i-1,j)=x_11(i-1,j)+1/2*minmod(x_11(i-1,j+1)-x_11(i-1,j),x_11(i-1,j)-x_11(i-1,j-1));
            hx12(i-1,j)=x_12(i-1,j)+1/2*minmod(x_12(i-1,j+1)-x_12(i-1,j),x_12(i-1,j)-x_12(i-1,j-1));
            hx13(i-1,j)=x_13(i-1,j)+1/2*minmod(x_13(i-1,j+1)-x_13(i-1,j),x_13(i-1,j)-x_13(i-1,j-1));
            hx22(i-1,j)=x_22(i-1,j)+1/2*minmod(x_22(i-1,j+1)-x_22(i-1,j),x_22(i-1,j)-x_22(i-1,j-1));
            hx23(i-1,j)=x_23(i-1,j)+1/2*minmod(x_23(i-1,j+1)-x_23(i-1,j),x_23(i-1,j)-x_23(i-1,j-1));
            hx33(i-1,j)=x_33(i-1,j)+1/2*minmod(x_33(i-1,j+1)-x_33(i-1,j),x_33(i-1,j)-x_33(i-1,j-1));
            hs11(i-1,j)=s_11(i-1,j)+1/2*minmod(s_11(i-1,j+1)-s_11(i-1,j),s_11(i-1,j)-s_11(i-1,j-1));
            hs12(i-1,j)=s_12(i-1,j)+1/2*minmod(s_12(i-1,j+1)-s_12(i-1,j),s_12(i-1,j)-s_12(i-1,j-1));
            hs13(i-1,j)=s_13(i-1,j)+1/2*minmod(s_13(i-1,j+1)-s_13(i-1,j),s_13(i-1,j)-s_13(i-1,j-1));
            hs22(i-1,j)=s_22(i-1,j)+1/2*minmod(s_22(i-1,j+1)-s_22(i-1,j),s_22(i-1,j)-s_22(i-1,j-1));
            hs23(i-1,j)=s_23(i-1,j)+1/2*minmod(s_23(i-1,j+1)-s_23(i-1,j),s_23(i-1,j)-s_23(i-1,j-1));
            hs33(i-1,j)=s_33(i-1,j)+1/2*minmod(s_33(i-1,j+1)-s_33(i-1,j),s_33(i-1,j)-s_33(i-1,j-1));
        end
        y1=y_1(i-1,j-1);
        y2=y_2(i-1,j-1);
        s11=s_11(i-1,j-1);
        s12=s_12(i-1,j-1);
        s13=s_13(i-1,j-1);
        s22=s_22(i-1,j-1);
        s23=s_23(i-1,j-1);
        s33=s_33(i-1,j-1);
        x11=x_11(i-1,j-1);
        x12=x_12(i-1,j-1);
        x13=x_13(i-1,j-1);
        x22=x_22(i-1,j-1);
        x23=x_23(i-1,j-1);
        x33=x_33(i-1,j-1);
        del1=delta1(i);del2=delta2(j);
        J=eval(Jacob_F_V); Der_1=eval(Der_F_del1); Der_2=eval(Der_F_del2);
        G=J\Der_1+J\Der_2;
        m_n(i,j)=m_n(i-1,j)+del_1*G(1)-del_1/del_2*(hm(i-1,j)-hm(i-1,j-1));
        y_1(i,j)=y_1(i-1,j)+del_1*G(2)-del_1/del_2*(hy1(i-1,j)-hy1(i-1,j-1));
        y_2(i,j)=y_2(i-1,j)+del_1*G(3)-del_1/del_2*(hy2(i-1,j)-hy2(i-1,j-1));
        x_11(i,j)=x_11(i-1,j)+del_1*G(4)-del_1/del_2*(hx11(i-1,j)-hx11(i-1,j-1));
        x_12(i,j)=x_12(i-1,j)+del_1*G(5)-del_1/del_2*(hx12(i-1,j)-hx12(i-1,j-1));
        x_13(i,j)=x_13(i-1,j)+del_1*G(6)-del_1/del_2*(hx13(i-1,j)-hx13(i-1,j-1));
        x_22(i,j)=x_22(i-1,j)+del_1*G(7)-del_1/del_2*(hx22(i-1,j)-hx22(i-1,j-1));
        x_23(i,j)=x_23(i-1,j)+del_1*G(8)-del_1/del_2*(hx23(i-1,j)-hx23(i-1,j-1));
        x_33(i,j)=x_33(i-1,j)+del_1*G(9)-del_1/del_2*(hx33(i-1,j)-hx33(i-1,j-1));
        s_11(i,j)=s_11(i-1,j)+del_1*G(10)-del_1/del_2*(hs11(i-1,j)-hs11(i-1,j-1));
        s_12(i,j)=s_12(i-1,j)+del_1*G(11)-del_1/del_2*(hs12(i-1,j)-hs12(i-1,j-1));
        s_13(i,j)=s_13(i-1,j)+del_1*G(12)-del_1/del_2*(hs13(i-1,j)-hs13(i-1,j-1));
        s_22(i,j)=s_22(i-1,j)+del_1*G(13)-del_1/del_2*(hs22(i-1,j)-hs22(i-1,j-1));
        s_23(i,j)=s_23(i-1,j)+del_1*G(14)-del_1/del_2*(hs23(i-1,j)-hs23(i-1,j-1));
        s_33(i,j)=s_33(i-1,j)+del_1*G(15)-del_1/del_2*(hs33(i-1,j)-hs33(i-1,j-1));
    end
end
%%
figure
surf(delta1,delta2,m_n)
    