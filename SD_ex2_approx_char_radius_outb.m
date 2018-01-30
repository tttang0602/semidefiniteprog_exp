load('boundary_rad1')
L = length(bound);
outb=zeros(L,2);
m=zeros(L,2);
h=0.001;
for i=1:L-1
    i
    k = bound(i,2)/bound(i,1);
    alp=1/sqrt(1+k^2);
    beta=abs(k)/sqrt(1+k^2);
    alp=alp*sign(bound(i,1));
    beta=beta*sign(bound(i,2));
    if bound(i,1)==0
        alp=0;
        beta=1;
    end
    [outb(i,1),outb(i,2),m(i,1),m(i,2)]=SD_ex2_approx_char_outerb(alp,beta,h);
end

% for i = 1:L-1
%     i
%     del1=bound_5(i,1);
%     del2=bound_5(i,2);
%     
%     while sum(abs(outb(i,:)))==0
%         if ne(bound_5(i,1),0)
%             k = bound_5(i,2)/bound_5(i,1);
%             alp=1/sqrt(1+k^2);
%             beta=abs(k)/sqrt(1+k^2);
%             alp=alp*sign(bound_5(i,1));
%             beta=beta*sign(bound_5(i,2));
%             del1=del1+del_v*sign(bound_5(i,1))*alp;
%             del2=del2+del_v*sign(bound_5(i,2))*beta;
%         else
%             del1=0;
%             del2=del2+del_v*sign(bound_5(i,2));
%         end
%         x=Exp2_solve(del1,del2);
%         if x(1)=='N'
%             i
%             outb(i,:)=[del1,del2];
%         elseif x<100
%             outb(i,:)=[del1,del2];
%         else
%             continue;
%         end
%     end
% end
%%
figure
plot(outb(:,1),outb(:,2))
       
