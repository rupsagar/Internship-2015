
clc
clear
global  U Baa k1 kk k2 Ba Jac b d d0 d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10

 Uval=50;
  T = 40000;
dt=0.1;
tau = 0:dt:T;
ics = [0 0 0 0 0 0 0 0];
% NT=length(tau);
l=1;
%%
%Taking Flutter values to compute [x](8x1)

  for n=1:length(Uval)
      U=Uval(n)
[t,x] = ode45('mylyapunov',[0:dt:T],ics);
z=eye(8);
expo_1=0; expo_2=0;expo_3=0;expo_4=0; expo_5 =0; expo_6=0; expo_7 =0; expo_8=0;
for i=500:size(x,1)
   
%Jacobian matrix (ref Lee 1998, AIAA Journal and Lee 1999 review paper)
% elements of jacobian B1...B8 and D1..D8 are defined below:
B1 = (d0*c6-c0*(d3+3*d4*transpose(x(i,1))*x(i,1)))/(d1*c0 -c1*d0);

B2 = (-c0*d2+d0*c3)/(d1*c0-c1*d0);

B3 = (d0*(c4+3*c5*transpose(x(i,3))*x(i,3))-c0*d6)/(d1*c0-c1*d0);

B4 = (d0*c2-c0*d5)/(d1*c0-c1*d0);

B5 = (d0*c7-c0*d7)/(d1*c0-c1*d0);

B6 = (d0*c8-c0*d8)/(d1*c0-c1*d0);

B7 = (d0*c9-c0*d9)/(d1*c0-c1*d0);

B8 = (d0*c10-c0*d10)/(d1*c0-c1*d0);

D1 = (-d1*c6+c1*(d3+3*d4*transpose(x(i,1))*x(i,1)))/(d1*c0-c1*d0);

D2 = (c1*d2-d1*c3)/(d1*c0-c1*d0);

D3 = (-d1*(c4+3*c5*transpose(x(i,3))*x(i,3))+c1*d6)/(d1*c0-c1*d0);

D4 = (-d1*c2+c1*d5)/(d1*c0-c1*d0);

D5 = (-d1*c7+c1*d7)/(d1*c0-c1*d0);

D6 = (-d1*c8+c1*d8)/(d1*c0-c1*d0);

D7 = (-d1*c9+c1*d9)/(d1*c0-c1*d0);

D8 = (-d1*c10+c1*d10)/(d1*c0-c1*d0);
   
    
    Jac = [0 1 0 0 0 0 0 0
       B1 B2 B3 B4 B5 B6 B7 B8
       0 0 0 1 0 0 0 0
       D1 D2 D3 D4 D5 D6 D7 D8
       1 0 0 0 -b 0 0 0
       1 0 0 0 0 -d 0 0
       0 0 1 0 0 0 -b 0
       0 0 1 0 0 0 0 -d]*dt;
X=x(i,:);
Y= transpose(X);
%% Gramshmidt orthonormalization to compute the 8 lyapunov exponents for a
% range of mean wind speed
    z=(eye(8)+Jac)*z;
    expo_1=expo_1+log(norm(z(:,1)));
    if i>1000
        exp_tim(l)=expo_1/(t(i)-t(500));l=l+1;
        i
    else
    end
    expo_2=expo_2+log(norm(z(:,2)));
    expo_3=expo_3+log(norm(z(:,3)));
    expo_4=expo_4+log(norm(z(:,4)));
    expo_5=expo_5+log(norm(z(:,5)));
    expo_6=expo_6+log(norm(z(:,6)));
    expo_7=expo_7+log(norm(z(:,7)));
    expo_8=expo_8+log(norm(z(:,8)));

        ortho_z=gramshmidt(z(:,1),z(:,2),z(:,3),z(:,4),z(:,5),z(:,6),z(:,7),z(:,8));
%     
        z=ortho_z;
end
expo_1=expo_1/(T-49.9);
expo_2=expo_2/(T-49.9);
expo_3=expo_3/(T-49.9);
expo_4=expo_4/(T-49.9);
expo_5=expo_5/(T-49.9);
expo_6=expo_6/(T-49.9);
expo_7=expo_7/(T-49.9);
expo_8=expo_8/(T-49.9);

Lyap(n,:)=[expo_1;expo_2;expo_3;expo_4;expo_5;expo_6;expo_7;expo_8]';
  end
 %plot(Uval,Lyap(:,1))
  