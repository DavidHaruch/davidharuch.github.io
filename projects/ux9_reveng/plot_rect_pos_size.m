function plot_rect_pos_size(x,y,xw,yw)

x1 = x;
y1 = y;
x2 = xw+x;
y2 = yw+y;

plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'k-')

end