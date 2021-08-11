function y1=rk4(f,t,x,v,h)
k=zeros(1,5);
k(1)=f(x,v,t);
k(2)=f(x+.5*h*v,v+.5*k(1)*h,t+.5*h);
k(3)=f(x+.5*h,y+.5*k(2)*h);
k(4)=f(x+h,y+h*k(3));
k(5)=h/6*(k(1)+2*k(2)+2*k(3)+k(4));
y1=y+k(5);