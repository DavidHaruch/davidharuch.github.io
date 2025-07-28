function plot_semicircle(x,y,d,theta1,theta2,style)

% theta = theta1:1:theta2;
theta = linspace(theta1,theta2,100);
x1 = x + (d/2).*sind(theta);
y1 = y + (d/2).*cosd(theta);
plot(x1,y1,style)
% plot(x,y,'*')

end