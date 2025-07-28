dia_range = linspace(0.3,0.4,50); % in mm

for z=1:length(dia_range)

% close all

% conductor_dia = 0.32; % in mm of only the conductor
conductor_dia = dia_range(z); % in mm of only the conductor

insulation_build = 0.01; % in mm of the insulation addition to diameter
wire_dia = conductor_dia + insulation_build;
wire_density = 2.7; % g/cc
wire_resistivity = 1.72e-8; % ohm meter at 20deg C
wire_temp_coeff = 3.93e-3; % % per deg
glue_density = 1.5; % g/cc
box_x = 4.5; % max width in mm
box_y = 11; % target height in mm
problem_length = 55; %mm


%% computation of the in-plane (wire stacking) measurements
% compute number of wires per first row
odd_row_qty = floor(box_x/wire_dia);
% compute number of wires per second row
offset_x = wire_dia/2;
offset_y = wire_dia*sind(60);
per_pair_height = wire_dia*(1+sind(60));

% figure(z)
% hold on
y_height_offset = 0;
j=1;
turns = 0;
full_row_number = 0;
full_1less_row_number = 0;
while y_height_offset+wire_dia/2 < box_y

    if mod(j,2)==1
        for i=1:odd_row_qty
%             plot_circle(wire_dia/2 + (i-1)*wire_dia,wire_dia/2+y_height_offset,wire_dia) % "odd row"
            turns = turns+1;
        end
        y_height_offset = y_height_offset + offset_y;
        full_row_number = full_row_number+1;
    else
        if 1
            for i=1:odd_row_qty-1
%             plot_circle(wire_dia/2 + (i-1)*wire_dia + offset_x,wire_dia/2+y_height_offset,wire_dia) % "even row"
            turns = turns+1;
            end
            full_1less_row_number = full_1less_row_number+1;
        else
            for i=1:odd_row_qty-1+1
%             plot_circle(wire_dia/2 + (i-1)*wire_dia + offset_x,wire_dia/2+y_height_offset,wire_dia) % "even row"
            turns = turns+1;
            end
            full_row_number = full_row_number+1;
        end
        y_height_offset = y_height_offset + offset_y;
    end
    j=j+1;
end

% plot_rect_pos_size(0,0,box_x,box_y)
% axis image
turns;
turns_arr(z) = turns;
actual_box_x = odd_row_qty*wire_dia;
actual_box_y = wire_dia+y_height_offset-offset_y;
% plot_rect_pos_size(0,0,box_x,actual_box_y)

% fill factor calculation
box_area = actual_box_x*actual_box_y; %mm^2
wire_area = turns*pi()*conductor_dia^2/4; % mm^2 area of conductor only
fill_factor(z) = wire_area/box_area;
% eff_density = fill_factor*wire_density + (1-fill_factor)*glue_density;

% total winding resistance calculation
per_side_resistance = turns*wire_resistivity*1000*(problem_length/(wire_area./turns)); % for one side of the coil in ohms
per_coil_resistance = 2*per_side_resistance; % ohms does not account for end turns

end_turn_avg_length = 28.7; % mm 
per_coil_resistance_w_ends = turns*wire_resistivity*1000*((2*problem_length+end_turn_avg_length*2)/(wire_area./turns));

Rph20(z) = per_coil_resistance_w_ends;
Rph25(z) = Rph20(z).*(1+wire_temp_coeff.*(25-20));

% % compute resistance, mass, and cost in common wire materials
% r_mat_arr = [1.72e-8 (2.6e-8*0.85 + 0.15*1.72e-8) 2.6e-8 1.59e-8]; % resistivity ohm*m in order copper, cca15, aluminum, silver
% rho_mat_arr = [8.96 (2.7*0.85 + 0.15*8.96) 2.7 10.49]; % density g/cc in order copper, cca15, aluminum, silver
% cost_mat_arr = [9.75 (2.67*0.85 + 0.15*9.75) 2.67 870]; % cost USD$/kg in order copper, cca15, aluminum, silver


% title({['turns = ' num2str(turns)];...
%     ['fill factor = ' num2str(fill_factor)];...
%     ['coil resistance  = ' num2str(per_coil_resistance_w_ends) 'Ω'];...
%     ['target X = ' num2str(box_x) 'mm ; ' 'actual X = ' num2str(actual_box_x) 'mm'];...
%     ['target Y = ' num2str(box_y) 'mm ; ' 'actual Y = ' num2str(actual_box_y) 'mm']...
%     })

end

figure(1)
subplot(1,2,1)
plot(dia_range,Rph20)
hold on
plot(dia_range,Rph25)

plot(dia_range,15.9.*ones(size(Rph25)),'k--')
xlabel('conductor diameter (mm)')
ylabel('Rph')
legend('@ 20deg C','@ 25deg C','target')


% AWG is crazy
% The n gauge wire diameter dn in millimeters (mm) is equal to 0.127mm times 92 raised to the power of 36 minus gauge number n, divided by 39.

dia_range_awg = mm2awg(dia_range);

subplot(1,2,2)
plot(dia_range_awg,Rph20)
hold on
plot(dia_range_awg,Rph25)

plot(dia_range_awg,15.9.*ones(size(Rph25)),'k--')
xlabel('conductor diameter (AWG)')
ylabel('Rph')
legend('@ 20deg C','@ 25deg C','target')
set(gca, 'XDir','reverse') % AWG goes from hi to low