% NASTOS VIKTOR 
%AEM 9297
tspan = [0 25];
y0 = [15 50];

b=1
[t,x]=ode45(@(t,x) func(t,x,b),tspan, y0);
plot(x(:,1), 'm')
title("X1 Time Response")
figure 
plot(x(:,2), 'm')
title("X2 Time Response")
figure
plot(x(:,1),x(:,2), 'k')
xlabel('X1-axis')
ylabel('X2-axis')
hold on

b=2
[t,x]=ode45(@(t,x) func(t,x,b),tspan, y0);
plot(x(:,1),x(:,2), 'c')

b=2.1
[t,x]=ode45(@(t,x) func(t,x,b),tspan,y0);
plot(x(:,1),x(:,2), 'y')
hold off

legend('b=1','b=2','b=2.1')



function xdot=func(t,x,b)
xdot(1)=x(2);
if(x(2)>=b)
     xdot(2)=(-((x(2)-b) - 2)^2 - 3 - 3*x(1))/3;
else 
     xdot(2)=(((x(2)-b)+2)^2 + 3 - 3*x(1))/3;
end

xdot=xdot';
end
