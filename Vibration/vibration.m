function vibration(problem_type,k_or_l,m,method,z1,x0,v0,tspan,f,n)
% VIBRATION plots the analytical and numerical solutions superimposed on
% one another for velocity and displacement of a free or forced vibration,
% in two separate figure windows. The input for 'Simple Pendulum' type must
% be in CGS unit. The angular displacement initial condition of 'Simple
% Pendulum' type must be in degree.

% problem_type='Spring Mass' or 'Simple Pendulum'
% k_or_l=spring stiffness in 'Spring Mass' or length in 'Simple Pendulum'
% m=mass of the body
% method='varying_initial','varying_none','varying_zeta' or 'varying_both'
% z1=vector containing the different zeta values
% x0=vector containing initial displacement condition
% v0=initial velocity condition
% tspan=time over which the numerical solution is calculated and the plots are done
% f=external periodic force
% n=factor defining the no. of times forcing frequency is of the natural frequency

% A choice of input can be "vibration('Spring Mass',1000,5,'varying_zeta',[0 .1 .25 .5 .75 1],.05,0,[0 5],0,0)"
switch problem_type 
    case 'Spring Mass'
        k=k_or_l;
        on=(k/m)^.5;
        disp_ylabel='Displacement(m)'; vel_ylabel='Velocity(m/s)';
    case 'Simple Pendulum'
        g=981; l=k_or_l; x0=x0*pi/180;
        on=(g/l)^.5;
        disp_ylabel={'Angular','Displacement(rad)'}; vel_ylabel={'Angular','Velocity(rad/s)'};
    otherwise
        error('Invalid Problem Type. The Valid Options Are "Spring Mass" & "Simple Pendulum"')
end
o=n*on;
n_z=numel(z1);

if f==0||n==0, heading1='Free'; else, heading1='Forced'; end

switch method
    case 'varying_initial'
    if n_z~=1||numel(x0)==1, error(['Method ',method,' takes one zeta and multiple initial conditions']); end
    case 'varying_none'
    if n_z~=1||numel(x0)~=1, error(['Method ',method,' takes only one zeta and initial condition']); end
    case 'varying_zeta'
    if n_z==1||numel(x0)~=1, error(['Method ',method,' takes multiple zeta and one initial condition']); end
    case 'varying_both'
    if n_z==1||numel(x0)==1, error(['Method ',method,' takes multiple zeta and initial conditions']); end
    otherwise
    error('Invalid Method');
end

for i=1:length(x0)    
    for j=1:length(z1)
        x0_1=x0(i);
        z=z1(j);
        ini=[x0_1 v0];
        c=2*m*on*z;
        switch problem_type 
            case 'Spring Mass'
                p=@(t,y)([y(2);-k/m*y(1)-c/m*y(2)+f/m*sin(o*t)]);
            case 'Simple Pendulum'
                z=z/l^2;
                p=@(t,y)([y(2);-g/l*y(1)-c/(m*l^2)*y(2)+f/(m*l^2)*sin(o*t)]);
        end
        a=-on*z;
        b=sqrt(1-z^2)*on;
        den=m*(on^2-o^2)^2+4*m*a^2*o^2;
        A=x0_1-2*f*a*o/den;      
        B=v0/b-x0_1*a/b+2*f*a^2*o/(b*den)-f*(on^2-o^2)*o/(b*den);
    
        [t,y]=ode45(p,tspan,ini);
    
        switch b
            case 0
                c1=x0_1-2*a*o*f/den;
                c2=-c1*a-f/den*(on^2-o^2)*o;
                f_dis=@(t)(exp(a*t)*(c1+c2*t)+f/den*((on^2-o^2)*sin(o*t)+2*a*o*cos(o*t)));
                f_vel=@(t)(exp(a*t)*(a*c1+c2*(1+a*t))+...
                    f/den*((on^2-o^2)*o*cos(o*t)-2*a*o^2*sin(o*t)));
            otherwise
                f_dis=@(t)(exp(a*t)*(A*cos(b*t)+B*sin(b*t))+...
                    f*((on^2-o^2)*sin(o*t)+2*a*o*cos(o*t))/den);
                f_vel=@(t)(exp(a*t)*((b*B+a*A)*cos(b*t)+(a*B-A*b)*sin(b*t))+...
                    f/den*((on^2-o^2)*o*cos(o*t)-2*a*o^2*sin(o*t)));
        end
        
        switch method 
            case {'varying_zeta','varying_both'}
                figure(2*i-1)
                subplot(ceil(n_z/2),2,j),fplot(f_dis,tspan,'r'),hold on
                subplot(ceil(n_z/2),2,j),plot(t,y(:,1))
                title(['\zeta=',num2str(z1(j))]);
                xlabel('Time(s)');ylabel(disp_ylabel);
                figure(2*i)
                subplot(ceil(n_z/2),2,j),fplot(f_vel,tspan,'r'),hold on
                subplot(ceil(n_z/2),2,j),plot(t,y(:,2))
                title(['\zeta=',num2str(z1(j))]);
                xlabel('Time(s)');ylabel(vel_ylabel);
            case {'varying_initial','varying_none'}
                figure(1),fplot(f_dis,tspan,'r'),hold on,plot(t,y(:,1))
                figure(2),fplot(f_vel,tspan,'r'),hold on,plot(t,y(:,2))
        end
    end
end
switch method
    case {'varying_zeta','varying_both'}
        for i=1:length(x0)
            figure(2*i-1),suptitle(['Displacement vs Time (x_0=',num2str(x0(i)),')'])
            legend('Analytical','Numerical');
            figure(2*i),suptitle(['Velocity vs Time (x_0=',num2str(x0(i)),')'])
            legend('Analytical','Numerical');
        end
    case {'varying_initial','varying_none'}
        if z1==0, heading2='Undamped'; else, heading2='Damped'; end
        figure(1)
        title({'Displacement vs Time',[heading1,' ',heading2,' Oscillation Of ',problem_type,' For \zeta=',num2str(z1)]});
        xlabel('Time(s)');ylabel(disp_ylabel);legend('Analytical','Numerical');
        figure(2)
        title({'Velocity vs Time',[heading1,' ',heading2,' Oscillation Of ',problem_type,' For \zeta=',num2str(z1)]});
        xlabel('Time(s)');ylabel(vel_ylabel);legend('Analytical','Numerical');
end