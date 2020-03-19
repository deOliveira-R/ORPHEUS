% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function specifies the geometry and nodalization of the fuel
% rod (cylindrical fuel column in cylindrical cladding surrounded 
% by coolant) as well as some other parameters of the unit cell 
% similar to the pressurized water reactor (PWR) core unit cell.

function input_and_initialize_PWR_like
  
  global g th fr
  

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input fuel rod geometry and nodalization --------------------------------
  g.nz = 10; % number of axial nodes

  g.fuel.rIn = 0; % inner fuel radius (m)
  g.fuel.rOut = 4.12e-3; % outer fuel radius (m)
  g.fuel.nr = 20;  % number of radial nodes in fuel

  g.clad.rIn = 4.22e-3; % inner clad radius (m)
  g.clad.rOut = 4.75e-3; % outer clad radius (m)
  g.clad.nr = 5;
  
  g.cool.pitch = 13.3e-3; % square unit cell pitch (m)
  g.cool.rOut = sqrt(g.cool.pitch^2/pi); % equivalent radius of the unit cell (m)

  g.dz0 = 0.3 * ones(g.nz,1); % height of the node (m)
  g.dzGasPlenum = 0.2; % height of the fuel rod gas plenum assuming it is empty (m)



% input average power rating in fuel --------------------------------------
  th.qLHGR0 = [0      10     1e20; ...        % time (s)
               200e2  200e2  200e2         ]; % linear heat generation rate (W/m)
    
% input fuel rod parameters -----------------------------------------------
  fr.clad.fastFlux = [0      10     1e20; ...        % time (s)
                      1e13   1e13   1e13          ]; % fast flux in cladding (1/cm2-s)

  fr.fuel.FGR = [0      10     1e20; ...        % time (s)
                 0.03   0.03   0.03          ]; % fission gas release (-)

  fr.ingas.Tplenum = 533; % fuel rod gas plenum temperature (K)
  fr.ingas.p0 = 1; % as-fabricated helium pressure inside fuel rod (MPa)
  fr.fuel.por = 0.05 * ones(g.nz,g.fuel.nr); % initial fuel porosity (-)

% input channel geometry --------------------------------------------------
  g.aFlow = 8.914e-5 * ones(g.nz,1); % flow area (m2)

% input channel parameters ------------------------------------------------
  th.mdot0_ = [0      10     1000; ...        % time (s)  
               0.3    0.3    0.3           ]; % flowrate (kg/s) 0.3
  th.p0 = 16; % coolant pressure (MPa)
  th.T0 = 533.0; % inlet temperature (K)



% INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize fuel geometry ------------------------------------------------
  g.fuel.dr0 = (g.fuel.rOut-g.fuel.rIn)/(g.fuel.nr-1); % fuel node radial thickness (m)
  g.fuel.r0 = (g.fuel.rIn : g.fuel.dr0 : g.fuel.rOut)'; % fuel node radius (m)
  g.fuel.r0_ = [g.fuel.rIn; interp1(g.fuel.r0,1.5:g.fuel.nr-0.5)'; g.fuel.rOut]; % fuel node boundary (m)
  g.fuel.a0_ = repmat(2*pi*g.fuel.r0_', g.nz, 1) .*  repmat(g.dz0, 1, g.fuel.nr+1); % XS area of fuel node boundary (m2)
  g.fuel.v0 = repmat(pi*diff(g.fuel.r0_.^2)', g.nz, 1) .*  repmat(g.dz0, 1, g.fuel.nr); % fuel node volume (m3)
  g.fuel.vFrac = (g.fuel.rOut^2 - g.fuel.rIn^2) / g.cool.rOut^2;

% initialize clad geometry ------------------------------------------------
  g.clad.dr0 = (g.clad.rOut-g.clad.rIn)/(g.clad.nr-1); % clad node radial thickness (m)
  g.clad.r0 = (g.clad.rIn : g.clad.dr0 : g.clad.rOut)'; % clad node radius (m)
  g.clad.r0_ = [g.clad.rIn; interp1(g.clad.r0,1.5:g.clad.nr-0.5)'; g.clad.rOut]; % clad node boundary (m)
  g.clad.a0_ = repmat(2*pi*g.clad.r0_', g.nz, 1) .*  repmat(g.dz0, 1, g.clad.nr+1); % XS area of clad node boundary (m2)
  g.clad.v0 = repmat(pi*diff(g.clad.r0_.^2)', g.nz, 1) .*  repmat(g.dz0, 1, g.clad.nr); % clad node volume (m3)
  g.clad.vFrac = (g.clad.rOut^2 - g.clad.rIn^2) ./ g.cool.rOut^2;

% initialize gap geometry -------------------------------------------------
  g.gap.dr0 = (g.clad.rIn - g.fuel.rOut)*ones(1:g.nz); % initial cold gap (m)
  g.gap.r0_ = (g.clad.rIn + g.fuel.rOut)/2; % average gap radius (m)
  g.gap.a0_ = 2*pi*g.gap.r0_*ones(g.nz, 1) .* g.dz0; % XS area of the mid-gap (m2)
  g.gap.vFrac = (g.clad.rIn^2 - g.fuel.rOut^2) ./ g.cool.rOut^2;
  
% initialize as-fabricated inner volumes and gas amount -------------------
  g.vGasPlenum = g.dzGasPlenum * pi*g.clad.rIn^2; % gas plenum volume (m3)
  g.vGasGap = g.dz0 .* pi*(g.clad.rIn^2 - g.fuel.rOut.^2); % gas gap volume (m3)
  g.vGasCentralVoid = g.dz0 .* pi*g.fuel.rIn.^2; % gas central void volume (m3)
  fr.ingas.muHe0 = fr.ingas.p0 * (g.vGasPlenum + sum(g.vGasGap + g.vGasCentralVoid)) / (8.31e-6 * 293); % as-fabricated gas amount inside fuel rod (mole)

% initialize gas gap status -----------------------------------------------
  g.gap.open = ones(g.nz,1);
  g.gap.clsd = zeros(g.nz,1);
  
% initialize fuel and clad total deformation components -------------------
  for i=1:3
      fr.fuel.eps0{i} = zeros(g.nz,g.fuel.nr);
      fr.clad.eps0{i} = zeros(g.nz,g.clad.nr);
  end

% initialize flow channel geometry ----------------------------------------
  g.volFlow = g.aFlow .* g.dz0; % volume of node (m3)
  g.areaHX = 2 * pi * g.clad.rOut * g.dz0; % heat exchange area(m2)(m2)
  g.dHyd = 4 * g.volFlow ./ g.areaHX; % hydraulic diameter (m)
  g.cool.vFrac = (g.cool.rOut^2 - g.clad.rOut^2) ./ g.cool.rOut^2;
  

% initialize thermal hydraulic parameters ---------------------------------

% Path to steam-water properties:
  path(path,['..' filesep '00.XSteam']);
  th.h0 = XSteam('h_pt', th.p0/10, th.T0-273)*1e3; % water enthalpy at core inlet (J/kg)
  th.h = th.h0*ones(g.nz,1); % initial enthalpy in nodes (kJ/kg)
  th.p = th.p0 * ones(g.nz,1); % initial pressure in nodes (MPa)     


end