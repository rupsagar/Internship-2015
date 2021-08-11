syms x r h
s=solve('h+r*x-x^3');
sdiff=diff(h+r*x-x^3);
figure
for i=-10:.05:10
    for j=1:size(s)
        s2=subs(s(j),{r,h},{i,3});       %h>0
        s3=double(s2);
        if isreal(s3)
            sdiff2=subs(sdiff,{x,r,h},{s2,i,3});   %h>0
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
title('Imperfect Bifurcation(h+rx-x^3=0)')