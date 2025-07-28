function plot_circle(x,y,d)

theta = 0:2:360;
x1 = x + (d/2).*sind(theta);
y1 = y + (d/2).*cosd(theta);
plot(x1,y1)
plot(x,y,'*')

end