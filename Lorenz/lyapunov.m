function lambda=lyapunov(T,dt,sigma,r,b)
% LYAPUNOV calculates all the Lyapunov exponents of Lorenz system.

% A choice of input can be "lyapunov(100,.01,10,28,8/3)"
lorenz=@(t,ini)([sigma*(ini(2)-ini(1));r*ini(1)-ini(2)-ini(1)*ini(3);...
    ini(1)*ini(2)-b*ini(3)]);
[t,y]=ode45(lorenz,0:dt:T,[1 2 3]);

lambda=zeros(1,3);
Q=eye(3);
l=zeros(1,3);
for i=1:numel(t)
    jcal=eye(3)+jacobian_lorenz(y(i,:))*dt;
    uqr=gram_schmidt(jcal*Q);
    Q=uqr{2};
    R=uqr{3};
    for k=1:3
        l(k)=log(abs((R(k,k))));
        lambda(k)=lambda(k)+l(k)/T;
    end
end
    function j=jacobian_lorenz(x)
%         JACOBIAN_LORENZ calculates the Jacobian of Lorenz system
        j=[-sigma,sigma,0;
            r-x(3),-1,-x(1);
            x(2),x(1),-b];
    end
    function uqr=gram_schmidt(a)
%         GRAM_SCHMIDT computes the QR decomposition of a matrix.
%         The output is a cell containing orthogonal(U), orthonormal(Q)
%         and upper triangular(R) matrices respectively.
        n=size(a,1);
        uu=zeros(n);
        qq=zeros(n);
        uu(:,1)=a(:,1);
        qq(:,1)=uu(:,1)/norm(uu(:,1));
        for ii=2:n
            uu(:,ii)=a(:,ii)-(a(:,ii)'*qq(:,1))*qq(:,1);
            for j=1:ii-1
                uu(:,ii)=uu(:,ii)-(uu(:,ii)'*qq(:,j))*qq(:,j);
            end
            qq(:,ii)=uu(:,ii)/norm(uu(:,ii));
        end
        rr=qq'*a;
        uqr={uu qq rr};
    end
end