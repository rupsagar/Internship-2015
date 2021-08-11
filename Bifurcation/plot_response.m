function plot_response(ini,r1)
% PLOT_RESPONSE plots the response of a bifurcation equation subjected to
% different initial conditions and different values of parameter r.

% 'saveas' command has been commented to prevent auto-saving of the figure

% A choice of input can be "plot_response([0 -1 -2],-1)"
tspan=[0 10];
col=['r' 'g' 'b' 'c' 'm' 'y' 'k'];

for i=1:numel(r1)
    figure
    for j=1:numel(ini)
        r=r1(i);
        p=@(t,x)(r+x^2);
        [t,y]=ode45(p,tspan,ini(j));
        plot(t,y,col(j),'linewidth',2),hold on
    end
    xlabel('T');ylabel('X');
    title({['X vs T response for R=',num2str(r)],'for Saddle Node Bifurcation (R+X^2)'})
    legend([repmat('X=',length(ini),1),num2str(ini(:))])
%     saveas(gcf,['Saddle_node',num2str(i),'.png']);
end