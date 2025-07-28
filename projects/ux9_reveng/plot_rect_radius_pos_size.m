function plot_rect_radius_pos_size(x,y,xw,yw,radius,color)

x1 = x;
y1 = y;
x2 = xw+x;
y2 = yw+y;

% plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'k-')

rectangle('Position',[x1 y1 xw yw],'Curvature',2.*radius./xw,'EdgeColor',color)

end