d0=-0.489;
d1=-0.49;
d2=-0.491;
d3=-0.492;
d4=-0.493;
y0=0.204305518472662;
y1=0.201107684592419;
y2=0.197482095669551;
y3=0.193115453265627;
y4=0.186900321003191;% 

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
% 


c0=-0.6;

for j=1:1000
    A=[1 (d0-c0)^(1/2) (d0-c0) (d0-c0)^(3/2);...
       1 (d1-c0)^(1/2) (d1-c0) (d1-c0)^(3/2);...
       1 (d2-c0)^(1/2) (d2-c0) (d2-c0)^(3/2);...
       1 (d3-c0)^(1/2) d3-c0 (d3-c0)^(3/2)];
   m=[y0 y1 y2 y3]';
   coeff=A\m;
   a0=coeff(1);
   a1=coeff(2);
   a2=coeff(3);
   a3=coeff(4);
    c0=c0-(a0+a1*(d4-c0)^(1/2)+a2*(d4-c0)+a3*(d4-c0)^(3/2)-y4)/(-a1/2*(d4-c0)^(-1/2)-a2-3/2*a3*(d4-c0)^(1/2));
    if c0>-0.49
        j
        break;
    end
end

a2=a2/2;
a3=a3/6;
