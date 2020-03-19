% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function initializes coolant parameters

function initializeCoolant

% Global structures
  global g clad cool

% channel flow area (m2)
  cool.aFlow0 = 8.914e-5;                                                  %INPUT
  cool.aFlow = cool.aFlow0 * ones(g.nz,1);

% volume of coolant node (m3)
  cool.volFlow = cool.aFlow .* g.dz0;

% heat exchange area(m2)
  cool.areaHX = 2 * pi * clad.rOut * g.dz0;

% hydraulic diameter (m)
  cool.dHyd = 4 * cool.volFlow ./ cool.areaHX;
  
% coolant inlet temperature (K) as a function of time (s)
  cool.T0t = [0.0    100.0  100.5  101.0  101.5  600.0; ...
              553.0  553.0  300.0  300.0  553.0  553.0];
  cool.T0 = cool.T0t(2,1);

% coolant pressure (MPa) as a function of time
  cool.p0t = [0.0      600.0; ...
              15.5     15.5];

% water enthalpy (J/kg) at core inlet
  cool.h0 = XSteam('h_pT', cool.p0t(2,1)*10, cool.T0-273.15)*1e3;

% initial coolant temperature (K), pressure (MPa) and enthalpy (J/kg)
  cool.T = cool.T0*ones(g.nz,1);
  cool.p = cool.p0t(2,1)*ones(g.nz,1);
  cool.h = cool.h0*ones(g.nz,1);

% coolant velocity at the channel inlet (m/s) as a function of time (s)
  cool.vel0t = [0.0    600.0; ...
                4.8    4.8];
end