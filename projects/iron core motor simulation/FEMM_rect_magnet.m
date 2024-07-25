function [] = FEMM_rect_magnet(x1,y1,x2,y2,material,dir,group)
%FEMM_rect draws a rectangle and adds desired properties
    mi_clearselected;
    mi_drawrectangle(x1,y1,x2,y2);
    dx=x2-x1;
    dy=y2-y1;
    mi_addblocklabel(x1+dx*0.05,y1+dy*0.05);
    mi_selectlabel(x1+dx*0.05,y1+dy*0.05);
    mi_setblockprop(material,0,0,0,dir,group,0);
    mi_clearselected;
end