% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function initializes neutronics parameters

function initializeNeutronics

% Global structures
  global kin fuel

% beta-effective
  kin.beff = [2.584E-4 1.520E-3 1.391E-3 3.070E-3 1.102E-3 2.584E-4 ];     %INPUT
  
% decay constants
  kin.lmb  = [0.013    0.032    0.119    0.318    1.403    3.929    ];     %INPUT

% prompt neutron lifetime
  kin.tLife = 20.e-6;                                                      %INPUT

% equilibrium concentrations of the delayed neutron precursors (1/m3)
  kin.cDNP = kin.beff * fuel.pow0 ./ (kin.tLife * kin.lmb);

% Doppler coefficient (1/K)
  kin.dopC = -2e-5;

% coolant temperature coefficient (1/K)
  kin.coolTemC = -20e-5;


end