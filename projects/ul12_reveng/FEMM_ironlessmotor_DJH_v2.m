%%FEMM 2 sided Voice Coil Motor
% D. Haruch May 2024

% openfemm(1);
openfemm;
newdocument(0);

problem_length = 55; %mm
mi_probdef(0,'millimeters','planar',1e-008,problem_length,30,0);
% mi_probdef(0,'millimeters','planar',1e-010,problem_length,45,0);

% materials
% magnet_grade = 'N35';
magnet_grade = 'N42';
% magnet_grade = 'N48';
% magnet_grade = 'N55';
back_iron_material = '1018 Steel';
wire_material = '1mm';
% wire_material = 'Copper';
halbach=0;
backironcut=0;

mi_getmaterial('Air');
mi_getmaterial(magnet_grade);
mi_getmaterial(back_iron_material);
mi_getmaterial(wire_material);
    
% define currents + and -

no_parallel_coils = 2; 
current = .5; % peak of sine, current in each coil if they were all in series

current = 35/2/.7071; % peak of sine, current in each coil if they were all in series

true_current = current*no_parallel_coils; % true current put into motor
current_rms = 0.7071*current; % rms of sine
true_current_rms = 0.7071*true_current; % rms of sine

mi_addcircprop('A',0,1);
mi_addcircprop('B',0,1);
mi_addcircprop('C',0,1);

mi_setcurrent('A',current);
mi_setcurrent('B',-current/2);
mi_setcurrent('C',-current/2);

%mi_addcircprop('A',current,0);
%%mi_addcircprop('B',-current/2,0);
%mi_addcircprop('C',-current/2,0);


% winding parameters
no_turns = 90;
phase_resistance = 0.64; %ohms
%phase_resistance = 4; %ohms

power = phase_resistance*(3/2)*true_current^2;

% all dimensions in mm unless otherwise stated

%coils

coil_thickness = 5; % along the y direction
coil_height = 10; % along the x direction
coil_gap = 7.25; % starts symetrically around 1/2 this value; gap between +ive and -ive bundles of turns
air_gap = 0.95+0.25; % per side distance from coil to magnet
coil_spacing = .75; % distance between coils

%boundary
boundary_x = 400;
boundary_y = 100;

%magnets
magnet_thickness = 4.5;
%magnet_thickness = 24;
magnet_height = 18;
% magnet_height = 15;

%back iron
backiron_thickness = 6;
backiron_height = 550/2; % 36.35
% backiron_height = 105;

%back iron cut
back_iron_cut_depth = 4; %mm
back_iron_cut_width = 12; %mm
back_iron_cut_angle = 0; %deg

% back_iron_cut_depth = 4.7; %mm
% back_iron_cut_width = 2; %mm
% back_iron_cut_angle = 72.5; %deg

%cyclic properties
magnetic_pitch = 21; %mm
magnetic_pitch = magnetic_pitch/2;
ns_unit_length = magnetic_pitch*2; %mm
coil_pitch = (coil_gap + coil_height*2 + coil_spacing);

%computed distances
coil_top_y = coil_thickness/2;
coil_bottom_y = -coil_top_y;
coil_middle_y = (coil_top_y + coil_bottom_y)/2;

start_magnet_y = coil_thickness/2 + air_gap;
end_magnet_y = start_magnet_y + magnet_thickness;
middle_magnet_y = (start_magnet_y + end_magnet_y)/2;
start_magnet_x = magnetic_pitch - magnet_height/2;
end_magnet_x = magnetic_pitch + magnet_height/2;
middle_magnet_x = (start_magnet_x+end_magnet_x)/2;
halbach_width = magnetic_pitch*2-magnet_height;
halbach_ratio = halbach_width./magnet_height;

start_backiron_y = start_magnet_y + magnet_thickness;
end_backiron_y = start_backiron_y + backiron_thickness;
middle_backiron_y = (start_backiron_y+end_backiron_y)/2;


%



%% draw points to pick air gap flux
 mi_addnode(backiron_height,0)
 mi_addnode(-backiron_height,0)

%% draw coils
phase_order = ['A','B','C'];
% for i=0:2
for i=0:11
    %this_phase = phase_order(i+1);
    this_phase = phase_order(mod(i,3)+1);
    offset = (i-1).*(coil_gap + coil_height*2 + coil_spacing);
    mi_drawrectangle(offset + coil_gap/2,coil_top_y,offset + coil_gap/2+coil_height,coil_bottom_y);
    mi_drawrectangle(offset + -coil_gap/2,coil_top_y,offset + -coil_gap/2-coil_height,coil_bottom_y);

    mi_addblocklabel((offset + coil_gap/2 + offset + coil_gap/2+coil_height)/2,coil_middle_y);
    mi_selectlabel((offset + coil_gap/2 + offset + -coil_gap/2-coil_height)/2,coil_middle_y);
%     mi_setblockprop(wire_material,0,0,strcat(this_phase,'+'),0,1,no_turns);
    mi_setblockprop(wire_material,0,0,strcat(this_phase),0,1,no_turns);
    mi_clearselected;

    mi_addblocklabel(offset - coil_gap/2 - coil_height/2,coil_middle_y);
    mi_selectlabel(offset - coil_gap/2 - coil_height/2,coil_middle_y);
%    mi_setblockprop(wire_material,0,0,strcat(this_phase,'-'),0,1,no_turns);
    mi_setblockprop(wire_material,0,0,strcat(this_phase),0,1,-no_turns);
    mi_clearselected;

    mi_selectrectangle(offset + coil_gap/2,coil_top_y,offset + coil_gap/2+coil_height,coil_bottom_y,4)
    mi_setnodeprop(this_phase,1)
    mi_selectrectangle(offset + -coil_gap/2,coil_top_y,offset + -coil_gap/2-coil_height,coil_bottom_y,4)
    mi_setnodeprop(this_phase,1) 
end

%% draw magnets
% loop to +x
% for i=-1:3
for i=-4:7
% for i=1:1
    offset = (i-1)*(magnetic_pitch*4);
    % top magnets
    mi_drawrectangle(start_magnet_x+offset,start_magnet_y,end_magnet_x+offset,end_magnet_y); % N magnet
    mi_drawrectangle(-start_magnet_x+offset,start_magnet_y,-end_magnet_x+offset,end_magnet_y); % S magnet
    mi_addblocklabel(middle_magnet_x+offset,middle_magnet_y);
    mi_selectlabel(middle_magnet_x+offset,middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,270,0,0);
    mi_clearselected;
    mi_addblocklabel(-middle_magnet_x+offset,middle_magnet_y);
    mi_selectlabel(-middle_magnet_x+offset,middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,90,0,0);
    mi_clearselected;

    % bottom magnets
    start_magnet_y = -start_magnet_y;
    end_magnet_y = -end_magnet_y;
    middle_magnet_y = -middle_magnet_y;

    mi_drawrectangle(start_magnet_x+offset,start_magnet_y,end_magnet_x+offset,end_magnet_y); % N magnet
    mi_drawrectangle(-start_magnet_x+offset,start_magnet_y,-end_magnet_x+offset,end_magnet_y); % S magnet
    mi_addblocklabel(middle_magnet_x+offset,middle_magnet_y);
    mi_selectlabel(middle_magnet_x+offset,middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,270,0,0);
    mi_clearselected;
    mi_addblocklabel(-middle_magnet_x+offset,middle_magnet_y);
    mi_selectlabel(-middle_magnet_x+offset,middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,90,0,0);
    mi_clearselected;


    if halbach
    if mod(i,2)==1
        halbachmod = 0;
%         halbachmod = 180;
    else
        halbachmod = 180
%         halbachmod = 0;
    end

    % top halbach magnet #2
    mi_drawrectangle(magnetic_pitch*2-start_magnet_x+offset,-start_magnet_y,magnetic_pitch*2+start_magnet_x+offset,-end_magnet_y); % S magnet
    mi_addblocklabel(magnetic_pitch*2+offset,-middle_magnet_y);
    mi_selectlabel(magnetic_pitch*2+offset,-middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,180+halbachmod,0,0);
    mi_clearselected;

    % top halbach magnet #1
    mi_drawrectangle(-start_magnet_x+offset,-start_magnet_y,start_magnet_x+offset,-end_magnet_y); % S magnet
    mi_addblocklabel(offset,-middle_magnet_y);
    mi_selectlabel(offset,-middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,0+halbachmod,0,0);
    mi_clearselected;

    % bottom halbach magnet #1
    mi_drawrectangle(-start_magnet_x+offset,start_magnet_y,start_magnet_x+offset,end_magnet_y); % S magnet
    mi_addblocklabel(offset,middle_magnet_y);
    mi_selectlabel(offset,middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,180+halbachmod,0,0);
    mi_clearselected;

    % bottom halbach magnet #2
    mi_drawrectangle(magnetic_pitch*2-start_magnet_x+offset,start_magnet_y,magnetic_pitch*2+start_magnet_x+offset,end_magnet_y); % S magnet
    mi_addblocklabel(magnetic_pitch*2+offset,middle_magnet_y);
    mi_selectlabel(magnetic_pitch*2+offset,middle_magnet_y);
    mi_setblockprop(magnet_grade,0,0,0,0+halbachmod,0,0);
    mi_clearselected;
    end
end

%% draw back irons
mi_drawrectangle(backiron_height,start_backiron_y,-backiron_height,end_backiron_y);
mi_addblocklabel(0,start_backiron_y*1.01);
mi_selectlabel(0,start_backiron_y*1.01);
mi_setblockprop(back_iron_material,0,0,0,0,0,0);
mi_clearselected;

if backironcut
% trapezoid back iron cut set to air
back_iron_cut_y = end_backiron_y - back_iron_cut_depth;
back_iron_cut_y_center = (back_iron_cut_y+end_backiron_y)/2;

add_width = tand(back_iron_cut_angle)*back_iron_cut_depth;
cut_half_width = back_iron_cut_width/2;
cut_bottom_width = back_iron_cut_width;
cut_top_width = cut_bottom_width + 2*add_width;
cut_top_half_width = cut_top_width/2;
area_of_trapezoid = back_iron_cut_depth*(cut_bottom_width + cut_top_width)*0.5;

% for m=-2:1
for m=-5:4
cx = (m*magnetic_pitch*2) + magnetic_pitch*1;
% cx = (m*magnetic_pitch*4) + magnetic_pitch*2;

mi_drawpolygon([cx-cut_top_half_width,end_backiron_y; ...
    cx+cut_top_half_width,end_backiron_y; ...
    cx+cut_half_width,back_iron_cut_y; ...
    cx-cut_half_width,back_iron_cut_y; ...
    cx-cut_top_half_width,end_backiron_y; ...
    ])
mi_drawpolygon([cx-cut_top_half_width,-end_backiron_y; ...
    cx+cut_top_half_width,-end_backiron_y; ...
    cx+cut_half_width,-back_iron_cut_y; ...
    cx-cut_half_width,-back_iron_cut_y; ...
    cx-cut_top_half_width,-end_backiron_y; ...
    ])


mi_addblocklabel(cx+1,back_iron_cut_y_center);
mi_selectlabel(cx+1,back_iron_cut_y_center);
mi_setblockprop('Air',0,0,0,0,0,0);
mi_clearselected;

mi_addblocklabel(cx+1,-back_iron_cut_y_center);
mi_selectlabel(cx+1,-back_iron_cut_y_center);
mi_setblockprop('Air',0,0,0,0,0,0);
mi_clearselected;

end
end

start_backiron_y = -start_backiron_y;
end_backiron_y = -end_backiron_y;
middle_backiron_y = -middle_backiron_y;

mi_drawrectangle(backiron_height,start_backiron_y,-backiron_height,end_backiron_y);
mi_addblocklabel(0,start_backiron_y*1.01);
mi_selectlabel(0,start_backiron_y*1.01);
mi_setblockprop(back_iron_material,0,0,0,0,0,0);
mi_clearselected;

%mi_saveas('abc123')
%% draw boundary

%draw boundary
mi_drawpolygon([boundary_x boundary_y; -boundary_x boundary_y; -boundary_x -boundary_y; boundary_x -boundary_y]);
mi_addblocklabel(boundary_x-1,boundary_y-1);
mi_selectlabel(boundary_x-1,boundary_y-1);
mi_setblockprop('Air',0,0,0,0,0,0);
mi_clearselected;
 
%mi_makeABC()
%%mi_addblocklabel(-50,-50);
%mi_selectlabel(-50,-50);
%%mi_setblockprop('Air',0,0,0,0,0,0);
%mi_clearselected;


mi_refreshview;
mi_zoomnatural;
mi_saveas(fullfile(pwd,'abc123.fem')); % this must be .fem

%% define commutation parameters
commutation_offset = -14;
commutate_flag = 0
%%
% pos_swp_test = -50:.1:50;
% 
% figure(99)
% for i=1:length(pos_swp_test)
%     pos = pos_swp_test(i);
%     elec_theta = 2*pi()*(pos-commutation_offset)/(magnetic_pitch*4);
%     elec_theta = elec_theta + (pi()/2); % 90deg off 
%     test_ia(i) = sin(elec_theta);
%     test_ib(i) = sin(elec_theta+(2*pi())/3);
%     test_ic(i) = sin(elec_theta+(4*pi())/3);
% end
% 
% hold off
% plot(pos_swp_test,test_ia)
% hold on
% plot(pos_swp_test,test_ib)
% plot(pos_swp_test,test_ic)

%% define range of positions swept

% pos_swp = linspace(-120,-120+42,20); % range of positions to sweep
%pos_swp = linspace(-102,-94,20); % range of positions to sweep
%pos_swp = -100; % range of positions to sweep
%%pos_swp = -97.5
% pos_swp=linspace(0,42,30);
pos_swp=29 - 42*3;
pos_swp=-76.7241;

% pos_swp=linspace(-97,-97+42,30);

force = [];
tic
for i=1:length(pos_swp)
    this_pos = pos_swp(i);

%         if (commutate_flag)
%     % compute new currents
%         elec_theta = 2*pi()*(this_pos-commutation_offset)/(magnetic_pitch*4);
%         elec_theta = elec_theta + (pi()/2); % 90deg off 
%         ia(i) = current*sin(elec_theta);
%         ib(i) = current*sin(elec_theta-(2*pi())/3);
%         ic(i) = current*sin(elec_theta-(4*pi())/3);
% 
% %         ia(i)=0
% %         ic(i) = 0;
% %         ib(i) = 0;
%         mi_setcurrent('A',ia(i));
%         mi_setcurrent('B',ib(i));
%         mi_setcurrent('C',ic(i));
%     end

    mi_seteditmode('group')
    mi_selectgroup(1)
    mi_movetranslate(this_pos,0)
    mi_analyse(1); % set param to 1 for no visible window
    mi_loadsolution();
    mo_groupselectblock(1)
    watts(i)=mo_blockintegral(6);
force(i)=mo_blockintegral(11); % lorentz force 
force_wst(i)=mo_blockintegral(18); % wieghted stress tensor
    mi_clearselected;
    mi_seteditmode('group')
    mi_selectgroup(1)
    mi_movetranslate(-this_pos,0)
%     mo_close()
end
toc
%%
kf = max(abs(force))/true_current; % N/Apk (KNS, parker?)
kf_rms = max(abs(force))/true_current_rms; % N/Arms (everyone else)
km = max(abs(force))/sqrt(power); % N/sqrt(W) (everyone else)
km_sq = km^2 % N^2/W (tecnotion)

%% grab the air gap flux density

ag_points = 5000;
ag_b = [];
x_ag_point = linspace(-backiron_height,backiron_height,ag_points);


for i=1:ag_points
 b1 = mo_getb(x_ag_point(i),0);
 ag_b(i) = b1(2);
end

figure(3)
hold on
plot(x_ag_point,ag_b)
xlabel('position (mm)')
ylabel('flux density B (T)')

%%
figure(2)
hold on
plot(pos_swp,force)
xlabel('Position (mm)')
ylabel('Force (N)')

figure(4)
hold on
plot(pos_swp,force./true_current_rms)
xlabel('Position (mm)')
ylabel('Force Constant (N/Arms)')

mo_close()