function [AWG] = mm2awg(mm)
%mm2awg converts dia in mm to AWG
AWG = 36 - 39.*(log(mm/0.127)/log(92));
end