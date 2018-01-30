clear sol
syms m y s11 s12 s22 x11 x12 x22

delta=-0.49:0.001:0;
num=length(delta);
%sol=zeros(num,1);
for n = 1:num
    
    del=delta(n);
    eqn1 = m-y == 0;
    eqn2 = -(2+del)*x11-(6+4*del)*x12+(1-3*del)*x22-1==0;
    eqn3 = 1+2*y-s11+del*(2+y)==0;
    eqn4 = -1+3*y-s12+del*(-1+2*y)==0;
    eqn5 = 2-y-s22+del*(3+3*y)==0;
    eqn6 = s11*x11+s12*x12==0;
    eqn7 = s11*x12+s12*x22==0;
    eqn8 = s12*x12+s22*x22==0;


    Sol =  solve(eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8);

    M = [Sol.m,Sol.y,Sol.s11,Sol.s12,Sol.s22,Sol.x11,Sol.x12,Sol.x22];
   
    N = size(Sol.m);
    k=1;
    for i = 1:N
        S=[Sol.s11(i),Sol.s12(i); Sol.s12(i),Sol.s22(i)];
        X=[Sol.x11(i),Sol.x12(i); Sol.x12(i),Sol.x22(i)];
        S_sdp=min(eig(S));
        X_sdp=min(eig(X));
        
        if (X_sdp>=0)&&(S_sdp>=0)&&(imag(M(i,1))==0)
            if k>1
                k 
            end    
            sol(n,:)=M(i,:);
            k=k+1;
        end

    end
end
sold=double(sol);
%save('Ex1_Sol_fine.mat','sold','delta')
%%
figure(1)
subplot(1,2,1)
plot(delta,real(sol),'LineWidth',2)
xlabel('\delta')
ylabel('Real(m)')
subplot(1,2,2)
plot(delta,imag(sol),'LineWidth',2)
xlabel('\delta')
ylabel('Imag(m)')
    
%%
clear all
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
ndel=51;
delta1=linspace(0,0.05,ndel);
%delta1=fliplr(linspace(-0.1,-0.05,ndel));
max_sol=zeros(ndel,ndel);
for n1=1:ndel
    del1=delta1(n1);
    n1 
    for n2=1:ndel
        n2
        del2=delta1(n2);
        
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

        Sol_2=solve(eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15);

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
                S_sdp=min(eig(real(S)));
                X_sdp=min(eig(real(X)));
                if ((X_sdp)>=-10^(-10))&&(S_sdp>=-10^(-10))
                    
                    sol_2(count)=M_2(i,1);
                    sol_full(count,:)=M_2(i,:);
                    count=count+1
                end
            end

        end
        
        max_sol(n1,n2)=max(sol_2);
        ind=find(sol_2==max_sol(n1,n2));
        sol(n1,n2,:)=sol_full(ind,:);
        %sol2=sol_2(find(sol_2(:,1)==max_sol),:);
    end
    
end
sol_del1=sol;
%%
save('Ex2_sol_001.mat','sol','delta1')
%%
figure(2)
mesh(delta1,delta1,max_sol)
%%
 load('Ex2_sol.mat')
 delta1=linspace(-0.05,0.05,21);
figure(4)
surf(delta1,delta1,sol(:,:,1))
xlabel('\delta_2')
ylabel('\delta_1')
  
    