% close all
conductor_dia = 0.64; % in mm of only the conductor
insulation_build = 0.065; % in mm of the insulation addition to diameter
wire_dia = conductor_dia + insulation_build;
wire_density = 8.96; % g/cc
wire_resistivity = 1.724102121803e-8; % ohm meter at 20deg C
wire_temp_coeff = 3.93e-3; % % per deg
glue_density = 1.5; % g/cc
box_x = 12.5; % UX coil % max width in mm 
box_y = 8; % UX coil % target height in mm

% box_x = 5.1; % max width in mm
% box_y = 10; % target height in mm

problem_length = 70; %mm
coil_thickness = 8; % along the y direction
coil_height = 13; % along the x direction
coil_gap = 11.5; % starts symetrically around 1/2 this value; gap between +ive and -ive bundles of turns
air_gap = 0.5; % per side distance from coil to magnet
coil_spacing = 5; % distance between coils

%% computation of the in-plane (wire stacking) measurements
% compute number of wires per first row
odd_row_qty = floor(box_x/wire_dia);
% compute number of wires per second row
offset_x = wire_dia/2;
offset_y = wire_dia*sind(60);
per_pair_height = wire_dia*(1+sind(60));

if odd_row_qty*wire_dia+offset_x>box_x
    % can't fit same # of odd row qty on teh even rows
    even_row_qty = odd_row_qty-1;
else
    % can fit same # of odd row qty on teh even rows
    even_row_qty = odd_row_qty;
end

figure(1)
subplot(1,3,1)
hold on
y_height_offset = 0;
j=1;
turns = 0;
odd_row_num = 0; % odd rows
even_row_num = 0; % even rows

while y_height_offset+wire_dia/2 < box_y

    if mod(j,2)==1 % odd rows
        for i=1:odd_row_qty
            plot_circle(wire_dia/2 + (i-1)*wire_dia,wire_dia/2+y_height_offset,wire_dia) % "odd row"
            turns = turns+1;
        end
        y_height_offset = y_height_offset + offset_y;
        odd_row_num = odd_row_num+1;
    else % even rows
        for i=1:even_row_qty
        plot_circle(wire_dia/2 + (i-1)*wire_dia + offset_x,wire_dia/2+y_height_offset,wire_dia) % "even row"
        turns = turns+1;
        end
        even_row_num = even_row_num+1;

        y_height_offset = y_height_offset + offset_y;
    end
    j=j+1;
end

plot_rect_pos_size(0,0,box_x,box_y)
axis image
if even_row_qty == odd_row_qty
    actual_box_x = odd_row_qty*wire_dia + offset_x;
else
    actual_box_x = odd_row_qty*wire_dia;
end
actual_box_y = wire_dia+y_height_offset-offset_y;
plot_rect_pos_size(0,0,box_x,actual_box_y)




%% plot the top down view (for odd layers)
figure(1)
subplot(1,3,2)
% axis image
hold on
title({'Top Down View of Coil';'ODD LAYERS'})

odd_layer_len = 0;
for j=1:odd_row_qty

    % plot solid lines for edges of wires
    plot_semicircle(0,problem_length,coil_gap + j*wire_dia*2,-90,90,'k-')
    plot_semicircle(0,0,coil_gap + j*wire_dia*2,-90,-270,'k-')
    plot([coil_gap/2 + j*wire_dia coil_gap/2 + j*wire_dia ],[0 problem_length],'k-')
    plot(-[coil_gap/2 + j*wire_dia coil_gap/2 + j*wire_dia ],[0 problem_length],'k-')

    
    this_end_turn_dia = -wire_dia + coil_gap + j*wire_dia*2;
    this_end_turn_length = pi()*this_end_turn_dia;
    
    odd_layer_len = odd_layer_len + this_end_turn_length; % accumulate the total end turn wire length per layer

    % plot dashed lines for center of wires
    plot_semicircle(0,problem_length,this_end_turn_dia,-90,90,'k--')
    plot_semicircle(0,0,this_end_turn_dia,-90,-270,'k--')
    plot([(coil_gap/2 + j*wire_dia -wire_dia/2) (coil_gap/2 + j*wire_dia -wire_dia/2) ],[0 problem_length],'k--')
    plot(-[(coil_gap/2 + j*wire_dia -wire_dia/2) (coil_gap/2 + j*wire_dia -wire_dia/2) ],[0 problem_length],'k--')
end

plot([coil_gap/2+box_x coil_gap/2+box_x],[0 problem_length],'r-')
plot([-coil_gap/2-box_x -coil_gap/2-box_x ],[0 problem_length],'r-')

plot([coil_gap/2 coil_gap/2],[0 problem_length],'r-')
plot([-coil_gap/2 -coil_gap/2],[0 problem_length],'r-')

plot_semicircle(0,problem_length,coil_gap + 2*box_x,-90,90,'r-')
plot_semicircle(0,0,coil_gap + 2*box_x,-90,-270,'r-')
plot_semicircle(0,problem_length,coil_gap,-90,90,'r-')
plot_semicircle(0,0,coil_gap,-90,-270,'r-')

%% plot the top down view (for even layers)
figure(1)
subplot(1,3,3)
% axis image
hold on
title({'Top Down View of Coil';'EVEN LAYERS'})

even_layer_len = 0;
for j=1:even_row_qty
% for j=1:1

    % plot solid lines for edges of wires
    plot_semicircle(0,problem_length,coil_gap + j*wire_dia*2 + offset_x*2,-90,90,'k-')
    plot_semicircle(0,0,coil_gap + j*wire_dia*2 + offset_x*2,-90,-270,'k-')
    plot_semicircle(0,problem_length,coil_gap + j*wire_dia*2 -wire_dia*2 + offset_x*2,-90,90,'k-')
    plot_semicircle(0,0,coil_gap + j*wire_dia*2 -wire_dia*2 + offset_x*2,-90,-270,'k-')

    plot([coil_gap/2 + j*wire_dia coil_gap/2 + j*wire_dia ],[0 problem_length],'k--')
    plot(-[coil_gap/2 + j*wire_dia coil_gap/2 + j*wire_dia ],[0 problem_length],'k--')
    
    plot([coil_gap/2 + j*wire_dia- wire_dia/2, coil_gap/2 + j*wire_dia - wire_dia/2],[0 problem_length],'k-')
    plot(-[coil_gap/2 + j*wire_dia+ wire_dia/2, coil_gap/2 + j*wire_dia + wire_dia/2],[0 problem_length],'k-')
    
    this_end_turn_dia = -wire_dia + coil_gap + j*wire_dia*2 + offset_x*2;
    this_end_turn_length = pi()*this_end_turn_dia; % full circle from both sides
    
    even_layer_len = even_layer_len + this_end_turn_length; % accumulate the total end turn wire length per layer

    % plot dashed lines for center of wires
    plot_semicircle(0,problem_length,this_end_turn_dia,-90,90,'k--')
    plot_semicircle(0,0,this_end_turn_dia,-90,-270,'k--')
    
end

plot([coil_gap/2+box_x coil_gap/2+box_x],[0 problem_length],'r-')
plot([-coil_gap/2-box_x -coil_gap/2-box_x ],[0 problem_length],'r-')
plot([coil_gap/2 coil_gap/2],[0 problem_length],'r-')
plot([-coil_gap/2 -coil_gap/2],[0 problem_length],'r-')
plot_semicircle(0,problem_length,coil_gap + 2*box_x,-90,90,'r-')
plot_semicircle(0,0,coil_gap + 2*box_x,-90,-270,'r-')
plot_semicircle(0,problem_length,coil_gap,-90,90,'r-')
plot_semicircle(0,0,coil_gap,-90,-270,'r-')

%% Total end turn length

odd_et_len = odd_layer_len*odd_row_num;
even_et_len = even_layer_len*even_row_num;
total_et_len = odd_et_len + even_et_len; % total length in mm of the end turns

%% Total non-end turn length

force_generating_wire_length = turns*2*problem_length; % straight wire segments, two sides
end_turn_ratio = total_et_len / (force_generating_wire_length + total_et_len); % lower is better
total_len = force_generating_wire_length + total_et_len; % total wire length in mm
total_len_m = total_len/1000; % length in m

%% fill factor calculation
box_area = actual_box_x*actual_box_y; %mm^2
wire_area = pi()*conductor_dia^2/4; % mm^2 area of conductor only
wire_area_m = wire_area*1e-6;
fill_factor = (wire_area*turns)/box_area;
eff_density = fill_factor*wire_density + (1-fill_factor)*glue_density;

%% Computation of the coil resistance

% = rho * L/A
per_coil_resistance_w_ends = wire_resistivity*(1000*total_len/wire_area);

Rph20 = per_coil_resistance_w_ends
Rph25 = Rph20*(1+wire_temp_coeff*(25-20))

Rph25/3

subplot(1,3,1)
title({['turns = ' num2str(turns)];...
    ['fill factor = ' num2str(fill_factor)];...
    ['coil R20  = ' num2str(Rph20) 'Ω'];...
    ['coil R25  = ' num2str(Rph25) 'Ω'];...
    ['target X = ' num2str(box_x) 'mm ; ' 'actual X = ' num2str(actual_box_x) 'mm'];...
    ['target Y = ' num2str(box_y) 'mm ; ' 'actual Y = ' num2str(actual_box_y) 'mm']...
    })
