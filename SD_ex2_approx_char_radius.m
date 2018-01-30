
del_v=0.001;

%==============================Positive delta1
dir = 1;
del = 0:del_v:2;
%=============================================

%=============================Negative delta1==
% dir = -1;
% del_v = dir*del_v;
% del = 0:del_v:-0.6;
%=============================================


L = length(del);
format long
load('Exp2_Jacob.mat')
load('Ex2_sol.mat') 
Ang = 180;
L1 = L; L2 = 1;
m_n = zeros(L1,L2);
y_1 = zeros(L1,L2);
y_2 = zeros(L1,L2);
x_11 = zeros(L1,L2);
x_12 = zeros(L1,L2);
x_13 = zeros(L1,L2);
x_22 = zeros(L1,L2);
x_23 = zeros(L1,L2);
x_33 = zeros(L1,L2);
s_11 = zeros(L1,L2);
s_12 = zeros(L1,L2);
s_13 = zeros(L1,L2);
s_22 = zeros(L1,L2);
s_23 = zeros(L1,L2);
s_33 = zeros(L1,L2);
Del1 = zeros(L1,1);
Del2 = zeros(L1,1);
Merr = zeros(Ang,1);
gap = 1;
bound = zeros(Ang/gap,2);
m_b = zeros(Ang/gap,1);

a_m=[];
a_y1=[];
a_y2=[];
a_s11=[];
a_s12=[];
a_s13=[];
a_s22=[];
a_s23=[];
a_s33=[];
a_x11=[];
a_x12=[];
a_x13=[];
a_x22=[];
a_x23=[];
a_x33=[];
a_del1=[];
a_del2=[];

ind1 = 11;
ind2 = 11;

m_n(1) = sol(ind1,ind2,1);
y_1(1) = sol(ind1,ind2,2);
y_2(1) = sol(ind1,ind2,3);
s_11(1) = sol(ind1,ind2,4);
s_12(1) = sol(ind1,ind2,5);
s_13(1) = sol(ind1,ind2,6);
s_22(1) = sol(ind1,ind2,7);
s_23(1) = sol(ind1,ind2,8);
s_33(1) = sol(ind1,ind2,9);
x_11(1) = sol(ind1,ind2,10);
x_12(1) = sol(ind1,ind2,11);
x_13(1) = sol(ind1,ind2,12);
x_22(1) = sol(ind1,ind2,13);
x_23(1) = sol(ind1,ind2,14);
x_33(1) = sol(ind1,ind2,15);

%===============Test all solution sheets======= 
% y0=real(y0);
% m_n(1) = y0(rt,1);
% y_1(1) = y0(rt,2);
% y_2(1) = y0(rt,3);
% s_11(1) = y0(rt,4);
% s_12(1) = y0(rt,5);
% s_13(1) = y0(rt,6);
% s_22(1) = y0(rt,7);
% s_23(1) = y0(rt,8);
% s_33(1) = y0(rt,9);
% x_11(1) = y0(rt,10);
% x_12(1) = y0(rt,11);
% x_13(1) = y0(rt,12);
% x_22(1) = y0(rt,13);
% x_23(1) = y0(rt,14);
% x_33(1) = y0(rt,15);
%==========
    
for i = 0:gap:360
i
    alp=cos(i/180*pi);
    beta=sin(i/180*pi);
    dlambdax=zeros(L,1);
    dlambdas=zeros(L,1);
    for k = 1:L-1
        
        Del1(k+1) = alp*del(k+1);
        Del2(k+1) = beta*del(k+1);
        clear m y1 y2 s1* s2* s3* x1* x2* x3*
        m = m_n(k);
        y1 = y_1(k);
        y2 = y_2(k);
        s11 = s_11(k);
        s12 = s_12(k);
        s13 = s_13(k);
        s22 = s_22(k);
        s23 = s_23(k);
        s33 = s_33(k);
        x11 = x_11(k);
        x12 = x_12(k);
        x13 = x_13(k);
        x22 = x_22(k);
        x23 = x_23(k);
        x33 = x_33(k);
        del1 = Del1(k);
        del2 = Del2(k);
        J = eval(Jacob_F_V);
        Der_1 = eval(Der_F_del1);
        Der_2 = eval(Der_F_del2);
        DerV_1=J\Der_1';
        DerV_2=J\Der_2';
        G = -alp*(DerV_1)-beta*(DerV_2);
        
        
        m_n(k+1) = m_n(k)+del_v*G(1);
        y_1(k+1) = y_1(k)+del_v*G(2);
        y_2(k+1) = y_2(k)+del_v*G(3);
        s_11(k+1) = s_11(k)+del_v*G(4);
        s_12(k+1) = s_12(k)+del_v*G(5);
        s_13(k+1) = s_13(k)+del_v*G(6);
        s_22(k+1) = s_22(k)+del_v*G(7);
        s_23(k+1) = s_23(k)+del_v*G(8);
        s_33(k+1) = s_33(k)+del_v*G(9);
        
        x_11(k+1) = x_11(k)+del_v*G(10);
        x_12(k+1) = x_12(k)+del_v*G(11);
        x_13(k+1) = x_13(k)+del_v*G(12);
        x_22(k+1) = x_22(k)+del_v*G(13);
        x_23(k+1) = x_23(k)+del_v*G(14);
        x_33(k+1) = x_33(k)+del_v*G(15);
        
        clear m y1 y2 s1* s2* s3* x1* x2* x3* index
        m = m_n(k+1);
        y1 = y_1(k+1);
        y2 = y_2(k+1);
        s11 = s_11(k+1);
        s12 = s_12(k+1);
        s13 = s_13(k+1);
        s22 = s_22(k+1);
        s23 = s_23(k+1);
        s33 = s_33(k+1);
        x11 = x_11(k+1);
        x12 = x_12(k+1);
        x13 = x_13(k+1);
        x22 = x_22(k+1);
        x23 = x_23(k+1);
        x33 = x_33(k+1);
        del1 = Del1(k+1);
        del2 = Del2(k+1);
        J = eval(Jacob_F_V);
        Der_1 = eval(Der_F_del1); Der_2 = eval(Der_F_del2);
        G1 = -alp*(J\Der_1')-beta*(J\Der_2');
        
        Err_0(k) = eval(f0);
        Err_1(k) = eval(f1);
        Err_2(k) = eval(f2);
        Err_3(k) = eval(f3);
        Err_4(k) = eval(f4);
        Err_5(k) = eval(f5);
        Err_6(k) = eval(f6);
        Err_7(k) = eval(f7);
        Err_8(k) = eval(f8);
        Err_9(k) = eval(f9);
        Err_10(k) = eval(f10);
        Err_11(k) = eval(f11);
        Err_12(k) = eval(f12);
        Err_13(k) = eval(f13);
        Err_14(k) = eval(f14);
        Err = [max(Err_0),max(Err_1),max(Err_2),max(Err_3),max(Err_4),...
            max(Err_5),max(Err_6),max(Err_7),max(Err_8),max(Err_9),...
            max(Err_10),max(Err_11),max(Err_12),max(Err_13),max(Err_14)];
        Merr((i+90)/gap) = max(Err);
        m_n(k+1) = m_n(k)+0.5*del_v*(G(1)+G1(1));
        y_1(k+1) = y_1(k)+0.5*del_v*(G(2)+G1(2));
        y_2(k+1) = y_2(k)+0.5*del_v*(G(3)+G1(3));
        s_11(k+1) = s_11(k)+0.5*del_v*(G(4)+G1(4));
        s_12(k+1) = s_12(k)+0.5*del_v*(G(5)+G1(5));
        s_13(k+1) = s_13(k)+0.5*del_v*(G(6)+G1(6));
        s_22(k+1) = s_22(k)+0.5*del_v*(G(7)+G1(7));
        s_23(k+1) = s_23(k)+0.5*del_v*(G(8)+G1(8));
        s_33(k+1) = s_33(k)+0.5*del_v*(G(9)+G1(9));
        x_11(k+1) = x_11(k)+0.5*del_v*(G(10)+G1(10));
        x_12(k+1) = x_12(k)+0.5*del_v*(G(11)+G1(11));
        x_13(k+1) = x_13(k)+0.5*del_v*(G(12)+G1(12));
        x_22(k+1) = x_22(k)+0.5*del_v*(G(13)+G1(13));
        x_23(k+1) = x_23(k)+0.5*del_v*(G(14)+G1(14));
        x_33(k+1) = x_33(k)+0.5*del_v*(G(15)+G1(15));
        
%===========================================================================
%                  Calculate derivative of lambda to delta

        m = m_n(k+1);
        y1 = y_1(k+1);
        y2 = y_2(k+1);
        s11 = s_11(k+1);
        s12 = s_12(k+1);
        s13 = s_13(k+1);
        s22 = s_22(k+1);
        s23 = s_23(k+1);
        s33 = s_33(k+1);
        x11 = x_11(k+1);
        x12 = x_12(k+1);
        x13 = x_13(k+1);
        x22 = x_22(k+1);
        x23 = x_23(k+1);
        x33 = x_33(k+1);
        del1 = Del1(k+1);
        del2 = Del2(k+1);
        J = eval(Jacob_F_V);
        Der_1 = eval(Der_F_del1); Der_2 = eval(Der_F_del2);
        DerV2_1=J\Der_1';
        DerV2_2=J\Der_2';
        G2 = -alp*(DerV2_1)-beta*(DerV2_2);
        
        X = double([x_11(k+1) x_12(k+1) x_13(k+1);...
            x_12(k+1) x_22(k+1) x_23(k+1);...
            x_13(k+1) x_23(k+1) x_33(k+1)]);
        S = double([s_11(k+1) s_12(k+1) s_13(k+1);...
            s_12(k+1) s_22(k+1) s_23(k+1);...
             s_13(k+1) s_23(k+1) s_33(k+1)]);
        if (max(imag(X(:)))<10^(-15)) || (max(imag(S(:)))<10^(-15))
            X=real(X);S=real(S);
            
        else
            index = k;
            break
        end
    
        [Vx,Dx,Wx]=eig(X);
        [Vs,Ds,Ws]=eig(S);
        [X_min(k),ox]=min(diag(Dx));
        [S_min(k),os]=min(diag(Ds));
        vx=Vx(:,ox);
        vs=Vs(:,os);
        wx=Wx(:,ox);
        ws=Ws(:,os);
        dx=[G2(10),G2(11),G2(12);...
            G2(11),G2(13),G2(14);...
            G2(12),G2(14),G2(15)];
        ds=[G2(4),G2(5),G2(6);...
            G2(5),G2(7),G2(8);...
            G2(6),G2(8),G2(9)];
        for p=1:3
            for q=1:3

                dlambdax(k)=dlambdax(k)+wx(p)*vx(q)/(vx'*wx)*dx(p,q);
                dlambdas(k)=dlambdas(k)+ws(p)*vs(q)/(vs'*ws)*ds(p,q);
            end
        end
%=====================Criteria Zero======Derivative of lambda==========        
        if abs(dlambdas(k))>1||abs(dlambdax(k))>1
            index=k
            break
        end
%======================================================================

%==========================Criterial One================SDP condition===========     
        
%         if ((X_min(k))<-10^(-5)) || ((S_min(k))<-10^(-5))
%             index = k
%             eig_v=[X_min S_min];
%             break
%         end
%===========================================================

% %======================Criterial Two===Derivative=======
%         if k>=2&&k<L-1
%             Der1_m=G2(1);
%             Der2_m=(DerV2_1(1)-DerV_1(1))/(alp*del_v)+(DerV2_2(1)-DerV_2(1))/(beta*del_v);
%             rat_m(k)=Der2_m/Der1_m;
%             Der1_y1=G2(2);
%                                                                                                                                                                                                                                                             
%             Der2_y1=(DerV2_1(2)-DerV_1(2))/(alp*del_v)+(DerV2_2(2)-DerV_2(2))/(beta*del_v);
%             rat_y1(k)=Der2_y1/(Der1_y1);
%             
%             Der1_y2=G2(3);
%             Der2_y2=(DerV2_1(3)-DerV_1(3))/(alp*del_v)+(DerV2_2(3)-DerV_2(3))/(beta*del_v);
%             rat_y2(k)=Der2_y2/(Der1_y2);
% 
%             Der1_s11=G2(4);
%             Der2_s11=(DerV2_1(4)-DerV_1(4))/(alp*del_v)+(DerV2_2(4)-DerV_2(4))/(beta*del_v);
%             rat_s11(k)=Der2_s11/(Der1_s11);
%             
%             Der1_s12=G2(5);
%             Der2_s12=(DerV2_1(5)-DerV_1(5))/(alp*del_v)+(DerV2_2(5)-DerV_2(5))/(beta*del_v);
%             rat_s12(k)=Der2_s12/(Der1_s12 );
%             
%             Der1_s13=G2(6);
%             Der2_s13=(DerV2_1(6)-DerV_1(6))/(alp*del_v)+(DerV2_2(6)-DerV_2(6))/(beta*del_v);
%             rat_s13(k)=Der2_s13/(Der1_s13 );
%             
%             Der1_s22=G2(7);
%             Der2_s22=(DerV2_1(7)-DerV_1(7))/(alp*del_v)+(DerV2_2(7)-DerV_2(7))/(beta*del_v);
%             rat_s22(k)=Der2_s22/(Der1_s22 );
%             Der1_s23=G2(8);
%             Der2_s23=(DerV2_1(8)-DerV_1(8))/(alp*del_v)+(DerV2_2(8)-DerV_2(8))/(beta*del_v);
%             rat_s23(k)=Der2_s23/(Der1_s23 );
%             Der1_s33=G2(9);
%             Der2_s33=(DerV2_1(9)-DerV_1(9))/(alp*del_v)+(DerV2_2(9)-DerV_2(9))/(beta*del_v);
%             rat_s33(k)=Der2_s33/(Der1_s33 );
%            
%             Der1_x11=G2(10);
%             Der2_x11=(DerV2_1(10)-DerV_1(10))/(alp*del_v)+(DerV2_2(10)-DerV_2(10))/(beta*del_v);
%             rat_x11(k)=Der2_x11/(Der1_x11 );
%             
%             Der1_x12=G2(11);
%             Der2_x12=(DerV2_1(11)-DerV_1(11))/(alp*del_v)+(DerV2_2(11)-DerV_2(11))/(beta*del_v);
%             rat_x12(k)=Der2_x12/(Der1_x12 );
%             
%             Der1_x13=G2(12);
%             Der2_x13=(DerV2_1(12)-DerV_1(12))/(alp*del_v)+(DerV2_2(12)-DerV_2(12))/(beta*del_v);
%             rat_x13(k)=Der2_x13/(Der1_x13 );
%             Der1_x22=G2(13);
%             Der2_x22=(DerV2_1(13)-DerV_1(13))/(alp*del_v)+(DerV2_2(13)-DerV_2(13))/(beta*del_v);
%             rat_x22(k)=Der2_x22/(Der1_x22 );
%             Der1_x23=G2(14);
%             Der2_x23=(DerV2_1(14)-DerV_1(14))/(alp*del_v)+(DerV2_2(14)-DerV_2(14))/(beta*del_v);
%             rat_x23(k)=Der2_x23/(Der1_x23 );
%             Der1_x33=G2(15);
%             Der2_x33=(DerV2_1(15)-DerV_1(15))/(alp*del_v)+(DerV2_2(15)-DerV_2(15))/(beta*del_v);
%             rat_x33(k)=Der2_x33/(Der1_x33 );
% 
% 
%             [rat(k),para]=max([rat_m(k),rat_y1(k),rat_y2(k),...
%                     rat_s11(k),rat_s12(k),rat_s13(k),rat_s22(k),rat_s23(k),rat_s33(k),...
%                     rat_x11(k),rat_x12(k),rat_x13(k),rat_x22(k),rat_x23(k),rat_x33(k)]);
%             if abs(rat(k))>=1057
%                 index = k
%                 break
%             end
%        end
%======================================================================
        
    end
    if (exist('index','var'))
        ind=index;
    else
        ind=L1-1;
    end
    m_b((i)/gap+1)=m_n(ind);
    bound((i)/gap+1,:) = [Del1(ind);Del2(ind)];
    %==========Save the approximation values =====================
%     a_m=[a_m;m_n(1:ind)];
%     a_y1=[a_y1;y_1(1:ind)];
%     a_y2=[a_y2;y_2(1:ind)];
%     a_s11=[a_s11;s_11(1:ind)];
%     a_s12=[a_s12;s_12(1:ind)];
%     a_s13=[a_s13;s_13(1:ind)];
%     a_s22=[a_s22;s_22(1:ind)];
%     a_s23=[a_s23;s_23(1:ind)];
%     a_s33=[a_s33;s_33(1:ind)];
%     a_x11=[a_x11;s_11(1:ind)];
%     a_x12=[a_x12;s_12(1:ind)];
%     a_x13=[a_x13;s_13(1:ind)];
%     a_x22=[a_x22;s_22(1:ind)];
%     a_x23=[a_x23;s_23(1:ind)];
%     a_x33=[a_x33;s_33(1:ind)];
%     a_del1=[a_del1;Del1(1:ind)];
%     a_del2=[a_del2;Del2(1:ind)];
   
    
    figure(9)
    plot3(Del1(1:ind),Del2(1:ind),m_n(1:ind))
    hold on
%     plot3(Del1(1:ind),Del2(1:ind),y_1(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),y_2(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),s_11(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),s_12(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),s_13(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),s_22(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),s_23(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),s_33(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),x_11(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),x_12(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),x_13(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),x_22(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),x_23(1:ind))
%     hold on
%     plot3(Del1(1:ind),Del2(1:ind),x_33(1:ind))
    %%===============Plot minimum eigenvalues================
%     figure(4)
%     plot3(Del1(1:ind),Del2(1:ind),X_min(1:ind))
%     hold on
%     xlabel('\delta_1')
%     zlabel('X_min')
%     %legend(['solution sheet ',num2str(rt-3)])
%     
%     figure(5)
%     plot3(Del1(1:ind),Del2(1:ind),S_min(1:ind))
%     hold on
%     xlabel('\delta_1')
%     zlabel('S_min')
%     %legend(['solution sheet ',num2str(rt-3)])

end


%save('Ex2_sol_app_p005.mat','a_*')
%save('boundary_rad1n','bound')
%%

xlabel('\delta_1')  
ylabel('\delta_2')
figure
plot(bound(1:end,1),bound(1:end,2),'LineWidth',3)
hold on
%%
plot(outb(1:72,1),outb(1:72,2),'LineWidth',3)