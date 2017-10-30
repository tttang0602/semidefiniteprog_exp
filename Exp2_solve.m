function [sol,All]=Exp2_solve(delta1,delta2)

%delta1=-0.045;delta2=0.045;
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



del1=delta1;

del2=delta2;
eqn1 = m-(y1-4*y2)==0;
eqn2 = trace((A1+del1*E11+del2*E21)*X_sym)-1==0;
eqn3 = trace((A2+del1*E12+del2*E22)*X_sym)+4==0;
 matr_stat = C-(y1*A1+y2*A2)+del1*(D1-y1*E11-y2*E12)+del2*(D2-y1*E21-y2*E22)-S_sym;
eqn4 = matr_stat(1,1)==0;
eqn5 = matr_stat(1,2)==0;
eqn6 = matr_stat(1,3)==0;
eqn7 = matr_stat(2,2)==0;
eqn8 = matr_stat(2,3)==0;
eqn9 = matr_stat(3,3)==0;

matr_comp=S_sym*X_sym;
eqn10 = matr_comp(1,1)==0;
eqn11 = matr_comp(1,2)==0;
eqn12 = matr_comp(1,3)==0;
eqn13 = matr_comp(2,2)==0;
eqn14 = matr_comp(2,3)==0;
eqn15 = matr_comp(3,3)==0;




eqns=[eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15];
vars=[m y1 y2 s11 s12 s13 s22 s23 s33 x11 x12 x13 x22 x23 x33];
Sol_2=solve(eqns,vars);

N = size(Sol_2.m);
M_2=double([Sol_2.m,Sol_2.y1,Sol_2.y2,...
    Sol_2.s11,Sol_2.s12,Sol_2.s13,Sol_2.s22,Sol_2.s23,Sol_2.s33,...
    Sol_2.x11,Sol_2.x12,Sol_2.x13,Sol_2.x22,Sol_2.x23,Sol_2.x33]);
count=1;

for i = 1:N
    if abs(imag(M_2(i,1)))<=10^(-15)
        
        S=double([Sol_2.s11(i),Sol_2.s12(i), Sol_2.s13(i);...
            Sol_2.s12(i),Sol_2.s22(i), Sol_2.s23(i);...
            Sol_2.s13(i),Sol_2.s23(i), Sol_2.s33(i)]);
        X=double([Sol_2.x11(i),Sol_2.x12(i), Sol_2.x13(i);...
            Sol_2.x12(i),Sol_2.x22(i), Sol_2.x23(i);...
            Sol_2.x13(i),Sol_2.x23(i), Sol_2.x33(i)]);
        S_sdp=min(eig(S));
        X_sdp=min(eig(X));
        if ((X_sdp)>=-10^(-10))&&(S_sdp>=-10^(-10))
            
            sol_2(count)=M_2(i,1);
            sol_full(count,:)=M_2(i,:);
            count=count+1;
        end
    end
    
end
if exist('sol_2','var')
    
    max_sol=max(sol_2);
    ind=find(sol_2==max_sol);
    sol=sol_full(ind,:);
else
    sol='No solution';
end
All=M_2;
% f0 = m-(y1-4*y2);
% f1 = trace((A1+del1*E11+del2*E21)*X_sym)-1;
% f2 = trace((A2+del1*E12+del2*E22)*X_sym)+4;
% f3 = matr_stat(1,1);
% f4 = matr_stat(1,2);
% f5 = matr_stat(1,3);
% f6 = matr_stat(2,2);
% f7 = matr_stat(2,3);
% f8 = matr_stat(3,3);
% 
% f9 = matr_comp(1,1);
% f10 = matr_comp(1,2);
% f11 = matr_comp(1,3);
% f12 = matr_comp(2,2);
% f13 = matr_comp(2,3);
% f14 = matr_comp(3,3);
% soltt=sol;
%     m=soltt(1);
%     y1=soltt(2);
%     y2=soltt(3);
%     s11=soltt(4);
%     s12=soltt(5);
%     s13=soltt(6);
%     s22=soltt(7);
%     s23=soltt(8);
%     s33=soltt(9);
%     x11=soltt(10);
%     x12=soltt(11);
%     x13=soltt(12);
%     x22=soltt(13);
%     x23=soltt(14);
%     x33=soltt(15);
% X=[x11 x12 x13; x12 x22 x23; x13 x23 x33];
% S=[s11 s12 s13 ; s12 s22 s23; s13 s23 s33];
% X_min=min(eig(X))
% S_min=min(eig(S))

%    errb=[eval(f0) eval(f1) eval(f2) eval(f3) eval(f4) eval(f5) eval(f6) eval(f7) eval(f8) eval(f9) eval(f10) eval(f11) eval(f12) eval(f13) eval(f14)]