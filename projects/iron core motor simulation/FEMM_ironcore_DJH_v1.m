%%FEMM 2 sided Voice Coil Motor
% D. Haruch May 2024

% openfemm(1);
openfemm;
newdocument(0);

problem_length = 55; %mm
mi_probdef(0,'millimeters','planar',1e-008,problem_length,30,0);

% materials
magnet_grade = 'N35';
back_iron_material = '1010 Steel';

% lam_material = 'M-19 Steel';
lam_material = 'Hiperco-50';
% lam_material = 'M-36 Steel';

coil_gap_material = 'Air';
halbach=0;
backironcut=0;

mi_getmaterial('Air');
mi_getmaterial(magnet_grade);
mi_getmaterial(back_iron_material);
mi_getmaterial('0.5mm');
mi_getmaterial(lam_material);

% define currents + and -
current = 10; % peak of sine
% current = 0; % peak of sine
current_rms = 0.7071*current; % rms of sine

% mi_addcircprop('A',current,1);
% mi_addcircprop('B',current,1);
mi_addcircprop('C',current,1);

% mi_addcircprop('A',current,1);
% mi_addcircprop('B',-current/2,1);
% mi_addcircprop('C',-current/2,1);

mi_addcircprop('A',0,1);
mi_addcircprop('B',0,1);
% mi_addcircprop('C',0,1);

no_turns = 1000;
phase_resistance = 16; %ohms
power = phase_resistance*(3/2)*current^2;

% all dimensions in mm unless otherwise stated

%% coils
coil_thickness = 20; % along the y direction
coil_height = 7; % along the x direction
coil_gap = 5; % starts symetrically around 1/2 this value; gap between +ive and -ive bundles of turns
coil_spacing = 1; % distance between coils
air_gap = 0.8; % distance from coil to magnet

%% cyclic properties
magnetic_pitch = 15; %mm north to south
coil_pitch = (coil_gap + coil_height*2 + coil_spacing);
no_magnets = 10; % even numbers only
no_3ph_coil_pairs = 3; % 1xx = 3 coils, 2x = 6 coils, etc...

%magnets
magnet_thickness = 4.5;
magnet_height = 13;% no more than magnetic pitch

%back iron
backiron_thickness = 3;
backiron_height = 160;
backiron_extra = (magnetic_pitch-magnet_height)/2; % extra backiron along travel dir in mm

%% draw coils and lamination stack

% draw phase A
offset = 0;
FEMM_rect_coil(offset+0,-air_gap,offset+coil_height,-air_gap-coil_thickness,'0.5mm','A',100,99)
FEMM_rect(offset+coil_height,-air_gap,offset+coil_gap+coil_height,-air_gap-coil_thickness,lam_material,99)
FEMM_rect_coil(offset+coil_gap+coil_height,-air_gap,offset+coil_gap+2*coil_height,-air_gap-coil_thickness,'0.5mm','A',-100,99)

FEMM_rect(offset+coil_gap+2*coil_height,-air_gap,offset+coil_spacing+coil_gap+2*coil_height,-air_gap-coil_thickness,coil_gap_material,99) % gap between coils


% draw phase B
offset = 1*(coil_spacing+coil_gap+2*coil_height);
FEMM_rect_coil(offset+0,-air_gap,offset+coil_height,-air_gap-coil_thickness,'0.5mm','B',100,99)
FEMM_rect(offset+coil_height,-air_gap,offset+coil_gap+coil_height,-air_gap-coil_thickness,lam_material,99)
FEMM_rect_coil(offset+coil_gap+coil_height,-air_gap,offset+coil_gap+2*coil_height,-air_gap-coil_thickness,'0.5mm','B',-100,99)

FEMM_rect(offset+coil_gap+2*coil_height,-air_gap,offset+coil_spacing+coil_gap+2*coil_height,-air_gap-coil_thickness,coil_gap_material,99) % gap between coils

% draw phase C
offset = 2*(coil_spacing+coil_gap+2*coil_height);
FEMM_rect_coil(offset+0,-air_gap,offset+coil_height,-air_gap-coil_thickness,'0.5mm','C',100,99)
FEMM_rect(offset+coil_height,-air_gap,offset+coil_gap+coil_height,-air_gap-coil_thickness,lam_material,99)
FEMM_rect_coil(offset+coil_gap+coil_height,-air_gap,offset+coil_gap+2*coil_height,-air_gap-coil_thickness,'0.5mm','C',-100,99)

% bottom of lam stack
FEMM_rect(0,-air_gap-coil_thickness,offset+coil_gap+2*coil_height,-air_gap-coil_thickness-10,lam_material,99)


%% draw back iron and magnets
FEMM_rect(-backiron_extra,magnet_thickness,no_magnets*magnetic_pitch-(magnetic_pitch-magnet_height)+backiron_extra,backiron_thickness+magnet_thickness,back_iron_material,1);

for i=0:no_magnets-1
    offset = i*magnetic_pitch;
    if mod(i,2)
        thisdir=90;
    else
        thisdir=270;
    end
    FEMM_rect_magnet(offset,0,magnet_height+offset,magnet_thickness,magnet_grade,thisdir,1);
end

%% draw boundary
% FEMM_rect(-50,-50,boundary_x,boundary_y,'Air',0);
mi_makeABC()
mi_addblocklabel(-50,-50);
mi_selectlabel(-50,-50);
mi_setblockprop('Air',0,0,0,0,0,0);

%%
%mi_saveas('abc123')
mi_refreshview;
mi_zoomnatural;
mi_saveas(fullfile(pwd,'abc123.fem')); % this must be .fem

%% define commutation parameters
commutation_offset = 35; % position in mm that maximizes force during DC load test (at this position, A=1, B=-0.5, C=-0.5)
% this is 90deg off from the zero force position found during typical
% commutation finding (woodpecker)
commutate_flag = 1; %do or not do sine commutation


%% sweep parameters

sweep_option = 4;
% 1 = position
% 2 = current
% 3 = cogging
% 4 = air gap

pos_swp = linspace(20,60,41); % in mm of stator movement
cur_swp = linspace(2,40,20); % in apk amps
airg_swp = linspace(0.5,15,15); % air gap in mm

%%
% pos_swp_test = 1:.1:100;
% figure(99)
% for i=1:length(pos_swp_test)
%     pos = pos_swp_test(i);
%     elec_theta = 2*pi()*(pos-commutation_offset)/(magnetic_pitch*2);
%     elec_theta = elec_theta + (pi()/2); % 90deg off 
% %     elec_theta = elec_theta - (pi()/2); % 90deg off 
%     test_ia(i) = sin(elec_theta);
% %     test_ib(i) = sin(elec_theta+(2*pi())/3);
% %     test_ic(i) = sin(elec_theta+(4*pi())/3);
%     test_ib(i) = sin(elec_theta+(4*pi())/3);
%     test_ic(i) = sin(elec_theta+(2*pi())/3);
% end
% 
% hold off
% plot(pos_swp_test,test_ia)
% hold on
% plot(pos_swp_test,test_ib)
% plot(pos_swp_test,test_ic)

%% define range of positions swept

if sweep_option == 1 || sweep_option ==  3
    fx = [];
    fz = [];
    ia = [];
    ib = [];
    ic = [];
    tic
    for i=1:length(pos_swp)
        this_pos = pos_swp(i);
        if (commutate_flag)
        % compute new currents
            elec_theta = 2*pi()*(this_pos-commutation_offset)/(magnetic_pitch*2);
            elec_theta = elec_theta + (pi()/2); % 90deg off
            if sweep_option == 3
                ia(i) = 0;
                ib(i) = 0;
                ic(i) = 0;
            else
                ia(i) = current*sin(elec_theta);
                ib(i) = current*sin(elec_theta+(4*pi())/3);
                ic(i) = current*sin(elec_theta+(2*pi())/3);
            end
            mi_setcurrent('A',ia(i));
            mi_setcurrent('B',ib(i));
            mi_setcurrent('C',ic(i));
        end
        % run simulation at new position
        mi_seteditmode('group')
        mi_selectgroup(99)
        mi_movetranslate(this_pos,0)
        mi_analyse(1); % set param to 1 for no visible window
        mi_resize(600,600)
        mi_loadsolution();
        mo_hidepoints()
        mo_resize(600,600)
        mo_showdensityplot(0,0,0,2,'mag')
        picture_name = strcat('images\',num2str(this_pos),'_ironcore.bmp')
        mo_savebitmap(picture_name)
        mo_groupselectblock(99)
        watts(i)=mo_blockintegral(6);
        fx(i)=mo_blockintegral(18);
        fz(i)=mo_blockintegral(19);
        mi_clearselected;
        mi_seteditmode('group')
        mi_selectgroup(99)
        mi_movetranslate(-this_pos,0)
    end
    toc
    %%
    % kf = max(abs(fx))/current; % N/Apk (KNS, parker?)
    % kf_rms = max(abs(fx))/current_rms; % N/Arms (everyone else)
    % km = max(abs(fx))/sqrt(power); % N/sqrt(W) (everyone else)
    % km_sq = km^2; % N^2/W (tecnotion)
    figure(1)
    subplot(1,2,1)
    hold on
    plot(pos_swp,fx)
    xlabel('Position (mm)')
    ylabel('Fx Force (N)')
    subplot(1,2,2)
    hold on
    plot(pos_swp,fz)
    xlabel('Position (mm)')
    ylabel('Fz Force (N)')
end

if sweep_option == 2
    fx = [];
    fz = [];
    ia = [];
    ib = [];
    ic = [];
    la = [];
    lb = [];
    lc = [];
    tic
    for i=1:length(cur_swp)
        this_cur = cur_swp(i);
%         this_pos = mean(pos_swp);
        this_pos = 50;
        if (commutate_flag)
        % compute new currents
            elec_theta = 2*pi()*(this_pos-commutation_offset)/(magnetic_pitch*2);
            elec_theta = elec_theta + (pi()/2); % 90deg off  
            ia(i) = this_cur*sin(elec_theta);
            ib(i) = this_cur*sin(elec_theta+(4*pi())/3);
            ic(i) = this_cur*sin(elec_theta+(2*pi())/3);
            mi_setcurrent('A',ia(i));
            mi_setcurrent('B',ib(i));
            mi_setcurrent('C',ic(i));
        end
        % run simulation at new position
        mi_seteditmode('group')
        mi_selectgroup(99)
        mi_movetranslate(this_pos,0)
        mi_analyse(1); % set param to 1 for no visible window
    %     mi_resize(600,600)
        mi_loadsolution();
    %     mo_hidepoints()
    %     mo_resize(600,600)
    %     mo_showdensityplot(0,0,0,2,'mag')
    %     picture_name = strcat('images\',num2str(this_pos),'_ironcore.bmp')
    %     mo_savebitmap(picture_name)
        mo_groupselectblock(99)
        watts(i)=mo_blockintegral(6);
        fx(i)=mo_blockintegral(18);
        fz(i)=mo_blockintegral(19);
        
        ca=mo_getcircuitproperties('A');
        cb=mo_getcircuitproperties('B');
        cc=mo_getcircuitproperties('C');

        la(i)=ca(3)/ca(1);
        lb(i)=cb(3)/cb(1);
        lc(i)=cc(3)/cc(1);

        mi_clearselected;
        mi_seteditmode('group')
        mi_selectgroup(99)
        mi_movetranslate(-this_pos,0)
    end
    toc
    %%
    % kf = max(abs(fx))/current; % N/Apk (KNS, parker?)
    % kf_rms = max(abs(fx))/current_rms; % N/Arms (everyone else)
    % km = max(abs(fx))/sqrt(power); % N/sqrt(W) (everyone else)
    % km_sq = km^2; % N^2/W (tecnotion)
    figure(2)
    subplot(1,4,1)
    hold on
    plot(cur_swp,fx)
    xlabel('Current (Apk)')
    ylabel('Fx Force (N)')
    subplot(1,4,2)
    hold on
    plot(cur_swp,fz)
    xlabel('Current (Apk)')
    ylabel('Fz Force (N)')

    subplot(1,4,3)
    hold on
    plot(cur_swp,1e3.*mean([la;lb;lc],1))
    xlabel('Current (Apk)')
    ylabel('Inductance (mH)')

    subplot(1,4,4)
    hold on
    plot(cur_swp(2:end),diff(fx)./diff(cur_swp))
    xlabel('Current (Apk)')
    ylabel('Force Constant (N/Apk)')

end


if sweep_option == 4
    pos_swp = 20;
    fx = [];
    fz = [];
    ia = [];
    ib = [];
    ic = [];
    la = [];
    lb = [];
    lc = [];
    tic
    for i=1:length(airg_swp)
        this_airg = airg_swp(i);
        this_airg_move = air_gap - this_airg
        if (commutate_flag)
        % compute new currents
            elec_theta = 2*pi()*(this_pos-commutation_offset)/(magnetic_pitch*2);
            elec_theta = elec_theta + (pi()/2); % 90deg off
            if sweep_option == 3
                ia(i) = 0;
                ib(i) = 0;
                ic(i) = 0;
            else
                ia(i) = current*sin(elec_theta);
                ib(i) = current*sin(elec_theta+(4*pi())/3);
                ic(i) = current*sin(elec_theta+(2*pi())/3);
            end
            mi_setcurrent('A',ia(i));
            mi_setcurrent('B',ib(i));
            mi_setcurrent('C',ic(i));
        end
        % run simulation at new position
        mi_seteditmode('group')
        mi_selectgroup(99)
        mi_movetranslate(pos_swp,this_airg_move)
        mi_analyse(1); % set param to 1 for no visible window
        mi_resize(600,600)
        mi_loadsolution();
        mo_hidepoints()
        mo_resize(600,600)
        mo_showdensityplot(0,0,0,2,'mag')
        picture_name = strcat('images_airgap\','img',num2str(i),'.bmp');
        mo_savebitmap(picture_name)
        mo_groupselectblock(99)
        watts(i)=mo_blockintegral(6);
        fx(i)=mo_blockintegral(18);
        fz(i)=mo_blockintegral(19);
        ca=mo_getcircuitproperties('A');
        cb=mo_getcircuitproperties('B');
        cc=mo_getcircuitproperties('C');
        la(i)=ca(3)/ca(1);
        lb(i)=cb(3)/cb(1);
        lc(i)=cc(3)/cc(1);
        mi_clearselected;
        mi_seteditmode('group')
        mi_selectgroup(99)
        mi_movetranslate(-pos_swp,-this_airg_move)
%         mi_movetranslate(0,-this_airg_move)
    end
    toc
    %%

    figure(1)
    subplot(1,3,1)
    hold on
    plot(airg_swp,fx./current)
    xlabel('Air Gap (mm)')
    ylabel('Force Constant (N/Apk)')
    subplot(1,3,2)
    hold on
    plot(airg_swp,fz)
    xlabel('Air Gap (mm)')
    ylabel('Fz Force (N)')

    subplot(1,3,3)
    hold on
    plot(airg_swp,1e3.*mean([la;lb;lc],1))
    xlabel('Air Gap (mm)')
    ylabel('Inductance (mH)')
end