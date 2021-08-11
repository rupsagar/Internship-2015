function poincare(eqn_fh,system,req_direction)
% POINCARE plots the Poincare Map for Rossler and Lorenz System. The
% function handle of the equation of the Poincare section, the name of the
% system and choice of direction are given as the input arguments.

% A choice of input can be "poincare(@(x,y,z)x-y,'lorenz','greater than zero')"
global x y z
syms x y z

eqn_str=char(eqn_fh); eqn_sym=sym(eqn_str(9:end));
diff_x=diff(eqn_sym,x);diff_y=diff(eqn_sym,y);diff_z=diff(eqn_sym,z);
grad=double([diff_x diff_y diff_z]);
dir_cos=grad/norm(grad); angle_rad=acos(dir_cos);
theta=acos(dir_cos(2)/sin(angle_rad(3)));

switch system
    case {'Rossler','rossler'}
        heading='Rossler'; tspan=[0 200]; a=.1; b=.1; c=14; ini1=[0 1 0];
        tangent_fh=@(t,ini)([(-ini(2)-ini(3));ini(1)+a*ini(2);b+ini(3).*(ini(1)-c)]);
        p=@(x1,ini)([1/(-ini(2)-ini(3));(x1+a*ini(2))/(-ini(2)-ini(3));...
            (b+ini(3)*(x1-c))/(-ini(2)-ini(3))]);
    case {'Lorenz','lorenz'}
        heading='Lorenz'; tspan=[0 50]; s1=10; r1=28; b1=8/3; ini1=[0 1 0];
        tangent_fh=@(t,ini)([s1*(ini(2)-ini(1));r1*ini(1)-ini(2)-ini(1)*ini(3);...
            ini(1)*ini(2)-b1*ini(3)]);
        p=@(x1,ini)([1/(s1*(ini(2)-x1));(r1*x1-ini(2)-x1*ini(3))/(s1*(ini(2)-x1));...
            (x1*ini(2)-b1*ini(3))/(s1*(ini(2)-x1))]);
    otherwise
        error('Invalid Choice Of System. The Valid Choices are: "Rossler" & "Lorenz"');
end
[t1,out1]=ode45(tangent_fh,tspan,ini1); t1_out1=[t1 out1];
eqn_val=eqn_fh(out1(:,1),out1(:,2),out1(:,3));

ii=0; lie=zeros(1000,4);
intersect=t1_out1(eqn_val==0,:);
if ~isempty(intersect)
    l1=check_direction(intersect(:,2:4));
    if any(l1), ii=ii+numel(l1==1); lie(1:ii,:)=intersect(l1,:); end
end

prod=eqn_val(1:end-1).*eqn_val(2:end);
cross1_pos2=[t1_out1(prod<0,:) find(prod<0)];
cross1=cross1_pos2(:,1:4); pos2=cross1_pos2(:,5);
lie2=zeros(size(cross1,1),4);
for i=1:size(cross1,1)
    index=pos2(i)+1;
    [x2,out2]=ode45(p,[cross1(i,2) out1(index,1)],[cross1(i,1) cross1(i,3) cross1(i,4)]);
    t_xyz=[out2(:,1) x2 out2(:,2) out2(:,3)];
    e=eqn_fh(t_xyz(:,2),t_xyz(:,3),t_xyz(:,4));
    k=find(abs(e)==min(abs(e)));
    lie2(i,:)=t_xyz(k(end),:);
    if check_direction(lie2(i,2:4))==1, ii=ii+1; lie(ii,:)=lie2(i,:); end
end

figure, hold on, xlabel('X'), ylabel('Y');
title({['Poincare Map For ',heading,' System'],['Through Poincare Section ',eqn_str(9:end),'=0'],['(direction ',req_direction,')']})
    if ~isempty(lie)
        lie(:,2:4)=my_rotate(lie(:,2:4),[3 1],[theta angle_rad(3)]);
        if dir_cos(2)<=0, lie(:,2:4)=my_rotate(lie(:,2:4),3,pi); else, view(0,-90); end
        plot(lie(:,2),lie(:,3),'b.');
    end
    if diff(minmax(lie(:,2)'))<.01, xlim([-50 50]); end
    if diff(minmax(lie(:,3)'))<.01, ylim([-50 50]); end

    function logical_dir=check_direction(xyz)
%         CHECK_DIRECTION checks the direction of the tangent to the phase space diagram at a particular point
        direction=sum(tangent_fh(0,xyz)'.*grad(1,:),2);
        switch req_direction
            case 'lesser than zero'
                logical_dir=direction<0;
            case 'greater than zero'
                logical_dir=direction>0;
            otherwise
                error('Invalid Choice Of Direction. The Valid Choices are: "lesser than zero" & "greater than zero"');
        end
    end
    function xyz=my_rotate(xyz,dir,angle)
%         MY_ROTATE rotates the plane to adjust the orientation 
        rx=@(phi)([1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)]);
        ry=@(phi)([cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)]);
        rz=@(phi)([cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1]);
        r={rx ry rz}; rot=eye(3);
        for m=1:numel(angle)
            rot=r{dir(m)}(angle(m))*rot;
        end
        for n=1:size(xyz,1)
            xyz(n,:)=(rot*xyz(n,:)')';
        end
    end
end