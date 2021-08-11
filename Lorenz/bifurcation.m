sigma=10;
r=28;
b=8/3;
rh=sigma*(sigma+b+3)/(sigma-b-1);
figure
% x0=sqrt(b*(rh-1));
% p=(x0)^2/(13.926-rh);
% q=(-x0)^2/(13.926-rh);


for r1=1:0.1:rh
x1=sqrt(b*(r1-1));
x2=-x1;
x8=0;
plot(r1,x1,'bo')
plot(r1,x2,'bo')
hold on
plot(r1,x8,'ro')
hold on
if r1>13.926
p=(x1)^2/(13.926-rh);
x3=sqrt(p*(r1-rh))+x1;
x4=-sqrt(p*(r1-rh))+x1;
q=(x2)^2/(13.926-rh);
x5=sqrt(q*(r1-rh))+x2;
x6=-sqrt(q*(r1-rh))+x2;
plot(r1,x3,'ro')
plot(r1,x4,'ro')
hold on
plot(r1,x5,'ro')
plot(r1,x6,'ro')
hold on
end
end
for r1=0:0.1:1
    x7=0;
    plot(r1,x7,'bo')
    hold on
end
for r1=rh:0.1:30
    x9=sqrt(b*(r1-1));
x10=-x9;
    x11=0;
    plot(r1,x9,'ro')
    plot(r1,x10,'ro')
    plot(r1,x11,'ro')
    hold on
end
legend('Stable','Unstable')
xlabel('r');ylabel('x');
title('Bifurcation Curve For Lorenz System','fontsize',16);
