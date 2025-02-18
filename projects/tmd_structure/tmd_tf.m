% tmd transfer function based code
% david john haruch 15 feb 2025

s = tf('s')

% main system
m = 10;
c = .1;
k = 100;

% TMD
mtmd = .5;
ctmd = .2;
ktmd = 5;

% nondimensionalized
mu = mtmd/m;
w1 = sqrt(k/m);
w2 = sqrt(ktmd/mtmd);
f = w1/w2;
zeta = ctmd/(2*sqrt(mtmd*ktmd));

zeta_opt = sqrt((3*mu)/(8*(1+mu)^3));
f_opt = 1/(1+mu);

w2_opt = w1*f_opt;
ktmd_opt = mtmd*w2_opt^2;
ctmd_opt = zeta_opt*2*sqrt(mtmd*ktmd_opt);

% solve for optimal tuning in un-nondimensionalized


G = 1/(m*s*s + c*s + k)
T = (mtmd*ctmd*s*s*s + mtmd*ktmd*s*s)/(mtmd*s*s + ctmd*s + ktmd)
Topt = (mtmd*ctmd_opt*s*s*s + mtmd*ktmd_opt*s*s)/(mtmd*s*s + ctmd_opt*s + ktmd_opt)
GTMD = (G)/(1+G*T)
GTMDopt = (G)/(1+G*Topt)

%close all
bodeplot(G)
hold on
bodeplot(GTMD)
bodeplot(GTMDopt)

legend('Baseline','TMD','Optimal TMD')