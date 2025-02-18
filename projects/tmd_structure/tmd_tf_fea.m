% tmd transfer function based code
% FEA mag/phase
% david john haruch 15 feb 2025
close all
s = tf('s');
p = bodeoptions('cstprefs');
p.FreqUnits = 'Hz';
p.MagUnits = 'abs';
p.MagScale = 'log';
%import FRF (units mm/N)
mag_csv = csvread("mag.csv");
phase_csv = csvread("phase.csv");
freq_csv = mag_csv(:,2);
mag = mag_csv(:,3)/1e3; % units N/m
phase = phase_csv(:,3);
plant_frd = [mag.* exp(1j*deg2rad(phase))];
Gfea = frd(plant_frd,freq_csv.*6.28);
G = tfest(Gfea,10,8);
% plot plant and estimate
figure(1)
bodeplot(Gfea,p)
hold on
bodeplot(G,p)
legend('FEA','10th order model fit')
% TMD parameters
mtmd = 0.2; %kg
ctmd = 50; %Ns/m
ktmd = 110e3; %N/m
mtmd2 = .05; %kg
ctmd2 = 50; %Ns/m
ktmd2 = 600e3; %N/m
mtmd3 = .05; %kg
ctmd3 = 5000; %Ns/m
ktmd3 = 1e3; %N/m
% transfer functions
T = (mtmd*ctmd*s*s*s + mtmd*ktmd*s*s)/(mtmd*s*s + ctmd*s + ktmd);
T2 = (mtmd2*ctmd2*s*s*s + mtmd2*ktmd2*s*s)/(mtmd2*s*s + ctmd2*s + ktmd2);
T3 = (mtmd3*ctmd3*s*s*s + mtmd3*ktmd3*s*s)/(mtmd3*s*s + ctmd3*s + ktmd3);
GTMD = (G)/(1+G*T);
GTMD2 = (G)/(1+G*(T+T2));
GTMD3 = (G)/(1+G*(T+T2+T3));

% plot results
figure(2)
bodeplot(G,p)
hold on
bodeplot(GTMD,p)
legend('Baseline','Plant + 1x TMD')

figure(3)
bodeplot(T,p)
hold on
bodeplot(T2,p)
bodeplot(T+T2,p)
bodeplot(T+T2+T3,p)

figure(4)
bodeplot(G,p)
hold on
bodeplot(GTMD,p)
bodeplot(GTMD2,p)
legend('Baseline','Plant + 1x TMD','Plant + 2x TMD')

figure(5)
bodeplot(G,p)
hold on
bodeplot(GTMD3,p)
legend('Baseline','Plant + 2x TMD + 1x LMD')