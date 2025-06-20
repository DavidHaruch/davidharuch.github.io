%%FEMM 2 sided Voice Coil Motor
% D. Haruch May 2024

% openfemm(1);
openfemm;
newdocument(0);

problem_length = 0.25*25.4; %
mi_probdef(0,'millimeters','planar',1e-008,problem_length,30,0);

% materials
% magnet_grade = 'N35';
magnet_grade = 'N42';
% magnet_grade = 'N48';
% magnet_grade = 'N55';
back_iron_material = '1018 Steel';

mi_getmaterial('Air');
mi_getmaterial(magnet_grade);
mi_getmaterial(back_iron_material);

%magnets

centermag_x = 0.0625*25.4;
centermag_y = 3.175*4;

sidemag_y = 2*0.125*25.4;
sidemag_x = 0.0625*25.4;
sidemag_gap_y = 0.5;
backiron_x = 2;

ag = 0.8;

sidemag_x_start = centermag_x/2 + ag;
sidemag_x_end = sidemag_x_start + sidemag_x;
sidemag_x_mid = (sidemag_x_start + sidemag_x_end)/2;
bi_x_end = sidemag_x_end+backiron_x;
bi_x_mid = (bi_x_end+sidemag_x_end)/2;


%%
mi_drawrectangle(-centermag_x/2,-centermag_y/2,centermag_x/2,centermag_y/2); % center magnet
mi_addblocklabel(0,0);
mi_selectlabel(0,0);
mi_setblockprop(magnet_grade,0,0,0,360,5,0); % group #5
mi_clearselected;
mi_selectrectangle(-centermag_x/2,-centermag_y/2,centermag_x/2,centermag_y/2)
mi_setnodeprop(magnet_grade,5)
mi_clearselected;

%%
mi_drawrectangle(sidemag_x_start,sidemag_gap_y,sidemag_x_end,sidemag_y); % R magnet 1
mi_addblocklabel(sidemag_x_mid,sidemag_y/2);
mi_selectlabel(sidemag_x_mid,sidemag_y/2);
mi_setblockprop(magnet_grade,0,0,0,360,1,0); % group #1
mi_clearselected;
mi_drawrectangle(sidemag_x_start,-sidemag_gap_y,sidemag_x_end,-sidemag_y); % R magnet 2
mi_addblocklabel(sidemag_x_mid,-sidemag_y/2);
mi_selectlabel(sidemag_x_mid,-sidemag_y/2);
mi_setblockprop(magnet_grade,0,0,0,180,1,0); % group #1
mi_clearselected;
mi_drawrectangle(sidemag_x_end,sidemag_y,bi_x_end,-sidemag_y); % backiron
mi_addblocklabel(bi_x_mid,0);
mi_selectlabel(bi_x_mid,0);
mi_setblockprop(back_iron_material,0,0,0,0,1,0); % group #1
mi_clearselected;


%%
sidemag_x_start = -sidemag_x_start;
sidemag_x_end = -sidemag_x_end;
sidemag_x_mid = -sidemag_x_mid;
bi_x_end = -bi_x_end;
bi_x_mid = -bi_x_mid;

mi_drawrectangle(sidemag_x_start,sidemag_gap_y,sidemag_x_end,sidemag_y); % L magnet 1
mi_addblocklabel(sidemag_x_mid,sidemag_y/2);
mi_selectlabel(sidemag_x_mid,sidemag_y/2);
mi_setblockprop(magnet_grade,0,0,0,360,1,0); % group #1
mi_clearselected;
mi_drawrectangle(sidemag_x_start,-sidemag_gap_y,sidemag_x_end,-sidemag_y); % L magnet 2
mi_addblocklabel(sidemag_x_mid,-sidemag_y/2);
mi_selectlabel(sidemag_x_mid,-sidemag_y/2);
mi_setblockprop(magnet_grade,0,0,0,180,1,0); % group #1
mi_clearselected;
mi_drawrectangle(sidemag_x_end,sidemag_y,bi_x_end,-sidemag_y); % backiron
mi_addblocklabel(bi_x_mid,0);
mi_selectlabel(bi_x_mid,0);
mi_setblockprop(back_iron_material,0,0,0,0,1,0); % group #1
mi_clearselected;
    
%% draw boundary


mi_makeABC()
mi_addblocklabel(3,8);
mi_selectlabel(3,8);
mi_setblockprop('Air',0,0,0,0,0,0);
mi_clearselected;


mi_refreshview;
mi_zoomnatural;
mi_saveas(fullfile(pwd,'abc123.fem')); % this must be .fem
%mi_analyse; % set param to 1 for no visible window
%mi_loadsolution();

%% define range of positions swept

pos_swp = linspace(5,-5,35); % range of positions to sweep
%pos_swp=-5
force_wst = [];
for i=1:length(pos_swp)
    this_pos = pos_swp(i);
    mi_seteditmode('group')
    mi_selectgroup(5)
    mi_movetranslate(0,this_pos)
    mi_analyse();
    mi_loadsolution();
    mo_groupselectblock(5)
    force_wst(i)=mo_blockintegral(19); % wieghted stress tensor
    force_wst_x(i)=mo_blockintegral(18); % wieghted stress tensor
    mo_clearblock;
    mi_deleteselected; 
    mi_clearselected;
    mi_seteditmode('group')
    mi_selectgroup(5)
    mi_movetranslate(0,-this_pos)
    mo_close()
end

force_wst

%% post processing +  stiffness calculation
air_gap=pos_swp;
air_gap_interp = linspace(min(air_gap),max(air_gap),1e3);
% fit curve to data
polyorder = 3;
p = polyfit(air_gap,force_wst,polyorder);
pder = polyder(p);
stiffness=polyval(pder,air_gap_interp); % dF/dz = stiffness (N/mm)

figure(1)
sgtitle({'Magnetic Gravity Compensator Linear Magnet Array'; ...
    'David J. Haruch';
    string(datetime)})
subplot(1,2,1)
hold on
plot(air_gap,force_wst)
xlabel('Z position (mm)')
ylabel('force (N)')
grid on
subplot(1,2,2)
hold on
plot(air_gap_interp, stiffness);
%plot(air_gap_interp,1.*ones(size(air_gap_interp)),'k-.')
xlabel('Z position  (mm)');
ylabel('Stiffness (N/mm)');
%legend('Magnetic Gravity Counterbalance','Elastic Spring/Flexure','Location','southeast')
grid on;

figure(2)
hold on
plot(air_gap,force_wst/9.81)
xlabel('Z position (mm)')
ylabel('weight (kg)')
grid on

figure(3)
hold on
plot(air_gap,force_wst)
plot(air_gap_interp,polyval(p,air_gap_interp),'k--')
legend('data','fit')

figure(4)
hold on
plot(air_gap,force_wst_x)