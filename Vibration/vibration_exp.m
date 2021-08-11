function vibration_exp(choice)
% VIBRATION_EXP plots the time response of experimental vibrations. The
% input is a vector containing choice of data to be plotted. The valid
% choices are, 0:Forced, 1:Free Damped, 2:Noise Added.

% A choice of input can be "vibration_exp([0 1 2])"
heading={'Forced','Free Damped','Noise Added'};
for i=1:length(choice)
    data=xlsread(['scope_',num2str(choice(i)),'.xls']);
    figure(choice(i)+1)
    plot(data(3:1002,1),data(3:1002,2)*9810/101.3)
    xlabel('Time(s)');ylabel('Acceleration(m/s^2)');
    title(['Experimental plot of Acceleration vs Time for ',heading{choice(i)+1},' Vibration'])
end