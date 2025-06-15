%%FEMM 2 sided Voice Coil Motor
% D. Haruch May 2024

% openfemm(1);
openfemm;
newdocument(0);
mi_probdef(0,'millimeters','axi',1e-008);

% materials
% magnet_grade = 'N35';
magnet_grade = 'N52';
% magnet_grade = 'N48';
% magnet_grade = 'N55';
back_iron_material = '1018 Steel';

mi_getmaterial('Air');
mi_getmaterial(magnet_grade);
mi_getmaterial(back_iron_material);
    

%boundary
boundary_r = 80;
boundary_h = 100;

%magnets
magnet_radius = 8*2;
magnet_height = 6*1.5;

ringmag_radius_ID = 6;
ringmag_radius_OD = 12*2;
ringmag_height = 4*1.5;
initgap = 2;

mi_drawrectangle(0,0,magnet_radius,-magnet_height); % N magnet
mi_addblocklabel(1,-1);
mi_selectlabel(1,-1);
mi_setblockprop(magnet_grade,0,0,0,270,5,0); % group #5
mi_clearselected;
mi_selectrectangle(0,0,magnet_radius,-magnet_height)
mi_setnodeprop(magnet_grade,5)
mi_clearselected;

mi_drawrectangle(0,1,magnet_radius+4,-magnet_height-4); % force region for N magnet weighted stress tensor
mi_addblocklabel(0.1,0.1);
mi_selectlabel(0.1,0.1);
mi_setblockprop('Air',0,0,0.5,0,5,0); % group #5
mi_clearselected;
mi_selectrectangle(0,1,magnet_radius+4,-magnet_height-4)
mi_setnodeprop('Air',5)
mi_clearselected;

mi_drawrectangle(ringmag_radius_ID,initgap,ringmag_radius_OD,initgap+ringmag_height); % S magnet
mi_addblocklabel(ringmag_radius_ID+1,initgap+1);
mi_selectlabel(ringmag_radius_ID+1,initgap+1);
mi_setblockprop(magnet_grade,0,0,0,90,0,0);
mi_clearselected;
    
%% draw boundary

%draw boundary
mi_drawpolygon([0 -boundary_h/2; 0 boundary_h/2; boundary_r boundary_h/2; boundary_r -boundary_h/2]);
mi_addblocklabel(boundary_r-1,0);
mi_selectlabel(boundary_r-1,0);
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
mi_analyse; % set param to 1 for no visible window
mi_loadsolution();

%% define range of positions swept

pos_swp = linspace(-3,-7,25); % range of positions to sweep
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
    mo_clearblock;
    mi_deleteselected; 
    mi_clearselected;
    mi_seteditmode('group')
    mi_selectgroup(5)
    mi_movetranslate(0,-this_pos)
    mo_close()
end

%% post processing +  stiffness calculation
air_gap=initgap-pos_swp;
air_gap_interp = linspace(min(air_gap),max(air_gap),1e3);
% fit curve to data
p = polyfit(air_gap,force_wst,3);
pder = polyder(p);
stiffness=polyval(pder,air_gap_interp); % dF/dz = stiffness (N/mm)

figure(1)
sgtitle({'Magnetic Gravity Compensator Ring+Round Magnet'; ...
    'David J. Haruch';
    string(datetime)})
subplot(1,2,1)
hold on
plot(initgap-pos_swp,force_wst)
xlabel('air gap (mm)')
ylabel('force (N)')
grid on
subplot(1,2,2)
hold on
plot(air_gap_interp, stiffness);
plot(air_gap_interp,1.*ones(size(air_gap_interp)),'k-.')
xlabel('air gap (mm)');
ylabel('Stiffness (N/mm)');
legend('Magnetic Gravity Counterbalance','Elastic Spring/Flexure','Location','southeast')
grid on;

figure(2)
hold on
plot(initgap-pos_swp,force_wst/9.81)
xlabel('air gap (mm)')
ylabel('weight (kg)')
grid on

figure(3)
hold on
plot(air_gap,force_wst)

plot(air_gap_interp,polyval(p,air_gap_interp),'k--')