syms x r 
s=solve('r*x-x^3');
sdiff=diff(r*x-x^3);
figure
for i=-10:0.05:10
    for j=1:size(s)
        s2=subs(s(j),i);
        s3=double(s2);
        if isreal(s3)
            sdiff2=subs(sdiff,{x,r},{s2,i});
            sdiff3=double(sdiff2);
            if sdiff3<0
                plot(i,s3,'r.')
                hold on
            elseif sdiff3>0
                plot(i,s3,'b.')
                hold on
%             else
%                 plot(i,s3,'c.')
%                 hold on
            end
        end
    end
end
%axis([-10 10 -5 5])
xlabel('r')
ylabel('x')
%legend('stable','unstable')
title('Supercritical Pitchfork Bifurcation(rx-x^3=0)')