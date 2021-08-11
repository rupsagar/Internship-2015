function lorenz(sigma,r,b,tspan,ini_condn,choice)
% LORENZ plots the X vs T, Y vs T, Z vs T, Z vs X, Z vs Y,
% 3D Phase Space Diagram and the Next Ampltude Map for the Lorenz system.
% There could be only one parameter having multiple entries.

% A choice of input can be "lorenz(10,28,8/3,[0 50],[0 1 0],[1 3])"

% 'saveas' commands have been commented to prevent auto-saving of figures
n_ini=size(ini_condn,1);
if n_ini==1, col='b'; else, col=['r' 'g' 'b' 'c' 'm' 'y']; end
leg=[repmat('[',n_ini,1),num2str(ini_condn),repmat(']',n_ini,1)];
min_max=zeros(1,n_ini*2);
for i=1:n_ini
    ini=ini_condn(i,:);
    for j=1:length(sigma)
        for k=1:length(r)
            for l=1:length(b)
                lorenz=@(t,ini)([sigma(j)*(ini(2)-ini(1));r(k)*ini(1)-ini(2)-ini(1)*ini(3);...
                    ini(1)*ini(2)-b(l)*ini(3)]);
                [t,out]=ode45(lorenz,tspan,ini);
                
                if any(choice==1)
                    if n_ini~=1, figure(1),hold on; else, figure();end
                    plot(t,out(:,1),col(i)), xlabel('T');ylabel('X');
                    title(['X vs T \sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))])
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'XT','.png'])
                end
                if any(choice==2)
                    if n_ini~=1, figure(2),hold on; else, figure();end
                    plot(t,out(:,2),col(i)), xlabel('T');ylabel('Y');
                    title(['Y vs T \sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))])
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'YT','.png'])
                end
                if any(choice==3)
                    if n_ini~=1, figure(3),hold on; else, figure();end
                    plot(t,out(:,3),col(i)), xlabel('T');ylabel('Z');
                    title(['Z vs T \sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))])
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'ZT','.png'])
                end
                if any(choice==4)
                    if n_ini~=1, figure(4),hold on; else, figure();end
                    plot(out(:,1),out(:,3),col(i)), xlabel('X');ylabel('Z');
                    title(['Z vs X \sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))])
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'ZX','.png'])
                end
                if any(choice==5)
                    if n_ini~=1, figure(5),hold on; else, figure();end
                    plot(out(:,2),out(:,3),col(i)), xlabel('Y');ylabel('Z');
                    title(['Z vs Y \sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))])
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'ZY','.png'])
                end
                if any(choice==6)
                    if n_ini~=1, figure(6),hold on; else, figure();end
                    plot3(out(:,1),out(:,2),out(:,3),col(i),'linewidth',2),xlabel('X');ylabel('Y');zlabel('Z'); grid on
                    title({'Phase Space Diagram for Lorenz System',['\sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))]})
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'3D_phase_space','.png'])
                end
                if any(choice==7)
                    pks=locate_peaks(out(:,3)');
                    if numel(pks)==0||numel(pks)==1
                        error(['Number of peaks:',num2str(numel(pks)),' for [',num2str(ini),']. Cannot find Next Amplitude Map']);
                    end
                    min_max(1,2*i-1:2*i)=minmax(pks);
                    if n_ini~=1, figure(7),hold on; else, figure();end
                    plot(pks(1:numel(pks)-1),pks(2:numel(pks)),[col(i),'.'])
                    if n_ini==1
                        if min_max(1)==min_max(2),min_max=[min_max(1)-10 min_max(2)+10];end
                        hold on,fplot('x',min_max,'k')
                        axis ([min_max min_max]),xlabel('Z_n');ylabel('Z_{n+1}');
                    end
                    title({'Next Amplitude Map',['Z_{n+1} vs Z_n \sigma=',num2str(sigma(j)),' r=',num2str(r(k)),' b=',num2str(b(l))]})
%                     saveas(gcf,['sigma',num2str(sigma(j)),'r',num2str(r(k)),'b',num2str(b(l)),'Zn1Zn','.png'])
                end
            end
        end
    end
end
if n_ini~=1
    for i=1:length(choice), figure(choice(i)),legend(leg);end
    if any(choice==7), figure(7)
        min_max_1=minmax(min_max);
        if min_max_1(1)==min_max_1(2),min_max_1=[min_max_1(1)-10 min_max_1(2)+10];end
        fplot('x',min_max_1,'k');
        axis ([min_max_1 min_max_1]),xlabel('Z_n');ylabel('Z_{n+1}');
    end
end

    function [pks,loc]=locate_peaks(x)
%         LOCATE_PEAKS calculates the local maximas in a given vector and also
%         finds out the positions of the maximas in the vector. The outputs are in
%         the form of row vectors.
    if ~any(size(x)==1), error('X should be a vector'); end
    iii=0;
    pks=zeros(size(x));
    loc=zeros(size(x));
    for ii=2:numel(x)-1
        if (x(ii)>x(ii+1))&&(x(ii)>x(ii-1))
            iii=iii+1;
            pks(iii)=x(ii);
            loc(iii)=ii;
        end
    end
    pks=pks(1:iii);
    loc=loc(1:iii);
    end
end