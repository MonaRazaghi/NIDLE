% Lower and upper bounds of 3219 reactions in 31 conditions
[lb,ub]=get_bounds();

load('iJO1366.mat')

% abundance matrix for 1026 reactions (catalyzed by homemeric enzymes) in 31 conditions
load('homomeric_abun.mat');

%solution includes the fluxes, y variable and indices of active reactions
[solution]=MILP(iJO1366,homomeric_abun,lb,ub);

%Calculate kapp based on flux and abundance
[Kapp,fluxes]=getkapp_MILP(homomeric_abun.reacind,homomeric_abun.abun,solution.flux,solution.y);


%Get the maximal apparent number
kmax=getkmax(Kapp,fluxes,homomeric_abun);