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
  cool.T0t = [0.0    200.0   210.0  5000.0; ...
              553.0  553.0   380.0  380.0];
  cool.T0 = cool.T0t(2,1);

% coolant pressure (MPa) as a function of time
  cool.p0t = [0.0      200.0    200.3    200.5    201.4    201.9    202.2    203.0   203.8   206.3   207.9   210.4   213.1   216.6   218.8   221.0   224.8   229.5   235.8   243.9   250.2   256.8   264.1   273.7   291.7   313.5   334.0   5000.0; ...
              15.50000 15.50000 14.65727 13.69425 12.27957 11.13598 10.41352 9.63108 8.54763 7.52432 7.52432 6.56115 5.35727 3.79223 2.70878 1.80575 1.20389 0.72230 0.54173 0.39122 0.36115 0.33108 0.39122 0.39122 0.33108 0.33108 0.33108 0.33108];

% water enthalpy (J/kg) at core inlet
  cool.h0 = XSteam('h_pT', cool.p0t(2,1)*10, cool.T0-273.15)*1e3;

% initial coolant temperature (K), pressure (MPa) and enthalpy (J/kg)
  cool.T = cool.T0*ones(g.nz,1);
  cool.p = cool.p0t(2,1)*ones(g.nz,1);
  cool.h = cool.h0*ones(g.nz,1);

% coolant velocity at the channel inlet (m/s) as a function of time (s)
  cool.vel0t = [0.0    200.0   201.0   202.    203.0   440.0    450.0   5000.0; ...
                4.8    4.8     2.0     1.0     0.002   0.002    0.1     0.1];
end