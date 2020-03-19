% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates properties for single-phase water, 
% single-phase steam or homogeneous steam-water mixture using 
% XSteam package


function [pro,Lsat,Vsat] = H2Oeq(p,h)


% convert enthalpy from J/kg to kJ/kg
  h = h/1e3;

% saturation temperature at p (K)
  pro.Tsat = XSteam('Tsat_p',p*10) + 273.15;

% saturation temperature   
  Lsat.T = pro.Tsat;
% liquid specific enthalpy at saturation (kJ/kg)
  Lsat.h = XSteam('hL_p',p*10);
% liquid specific heat at saturation (kJ/kg-K)
  Lsat.c_p = XSteam('CpL_p',p*10);
% liquid density at saturation (kg/m3)
  Lsat.rho = XSteam('rhoL_p',p*10);
% liquid surface tension at saturation (Pa)
  Lsat.sig = XSteam('st_p',p*10);
% liquid thermal conductivity at saturation (W/m-K)
  Lsat.k = XSteam('tcL_p',p*10);
% liquid dynamic viscosity at saturation (Pa-s)
  Lsat.mu = XSteam('my_pT',p*10,(pro.Tsat-273.15)-1);
% liquid kinematic viscosity at saturation (m2/s)
  Lsat.nu = Lsat.mu / Lsat.rho;

% steam enthalpy at saturation (kJ/kg)
  Vsat.h = XSteam('hV_p',p*10);
% steam specific isobaric heat at saturation (kJ/kg-K)
  Vsat.c_p = XSteam('CpV_p',p*10);
% steam density (kg/m3)
  Vsat.rho = XSteam('rhoV_p',p*10);
% steam thermal conductivity at saturation (W/m-K)
  Vsat.k = XSteam('tcV_p',p*10);
% steam dynamic viscosity at saturation (Pa-s)
  Vsat.mu = XSteam('my_pT',p*10,(pro.Tsat-273.15)+1);
% steam kinematic viscosity at saturation (m2/s)
  Vsat.nu = Vsat.mu / Vsat.rho;

% mixture quality (-)
  pro.x = (h - Lsat.h) / (Vsat.h - Lsat.h);

 % single phase liquid
  if pro.x < 0
   % temperature (K)
	 pro.T = XSteam('T_ph',p*10,h)+273.15;
   % void (-)
     pro.void = 0;
   % liquid density (kg/m3)
     pro.rho = XSteam('rho_ph',p*10,h);
   % liquid thermal conductivity (W/m-K)
     pro.k = XSteam('tc_ph',p*10,h);
   % liquid dynamic viscosity (Pa-s)
     pro.mu = XSteam('my_ph',p*10,h);
   % liquid kinematic viscosity (m2/s)
     pro.nu = pro.mu / pro.rho;
   % liquid specific isobaric heat (kJ/kg-K)
     pro.c_p = XSteam('Cp_ph',p*10,h);

% single phase steam
  elseif pro.x >= 1
   % steam temperature(K)
	 pro.T = XSteam('T_ph',p*10,h)+273.15;
   % void (-)
     pro.void = 1; 
   % steam density (kg/m3)
     pro.rho = XSteam('rho_ph',p*10,h);
   % steam thermal conductivity (W/m-K)
     pro.k = XSteam('tc_ph',p*10,h);
   % steam dynamic viscosity (Pa-s)
     pro.mu = XSteam('my_ph',p*10,h); 
   % steam kinematic viscosity (m2/s)
     pro.nu = pro.mu / pro.rho;
   % steam specific isobaric heat (kJ/kg-K)
     pro.c_p = XSteam('Cp_ph',p*10,h);

% two phase mixture
  else
   % steam-water mixture temperature(K)
     pro.T = pro.Tsat;
   % steam-water mixture void fraction (-)
     pro.void = XSteam('vx_ph',p*10,h);
   % steam-water mixture density at saturation (kg/m3)
     pro.rho = Lsat.rho * (1 - pro.void) + Vsat.rho * pro.void;
   % steam-water mixture thermal conductivity at saturation (W/m-K)
     pro.k = Lsat.k * (1 - pro.x) + Vsat.k * pro.x;
   % steam-water mixture dynamic viscosity at saturation (Pa-s)
     pro.mu = Lsat.mu * (1 - pro.x) + Vsat.mu * pro.x;
   % steam-water kinematic viscosity at tk (m2/s)
     pro.nu = pro.mu / pro.rho;
   % steam-water mixture specific isobaric heat (kJ/kg-K)
     pro.c_p = Lsat.c_p * (1 - pro.x) + Vsat.c_p * pro.x;
  end
  pro.p = p;
% convert specific enthalpy from kJ/kg to J/kg
  pro.h = h * 1e3;
% convert specific heat from kJ/kg-K to J/kg-K
  pro.c_p = pro.c_p * 1e3;
% convert liquid enthalpy at saturation from kJ/kg to J/kg
  Lsat.h = Lsat.h * 1e3;
% convert steam enthalpy at saturation from kJ/kg to J/kg
  Vsat.h = Vsat.h * 1e3;
end