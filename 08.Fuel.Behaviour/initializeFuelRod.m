% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function initializes fuel rod parameters

function initializeFuelRod

  clear all;

% Global structures
  global g fuel clad gap ingas

%--------------------------------------------------------------------------
% FUEL ROD MATERIAL PROPERTIES
%--------------------------------------------------------------------------
% correlations for material properties (in folder 00.Lib)
  [fuel, gap, clad] = matpro(fuel, gap, clad);

%-------------------------------------------------------------------------- 
% Initialize fuel

% initial height of the fuel (m)
  g.h0 = 3;                                                                %INPUT

% inner fuel radius (m) as fabricated
  fuel.rIn = 0;                                                            %INPUT

% outer fuel radius (m) as fabricated
  fuel.rOut = 4.12e-3;                                                     %INPUT

% number of radial nodes in fuel
  fuel.nr = 30;                                                            %INPUT

% initial porosity (-) assumed constant
  fuel.por = 0.05;                                                         %INPUT
  
% fission gas release (-) assumed constant
  fuel.FGR = 0.06;                                                         %INPUT

% linear heat generation rate (W/m) assumed constant
  fuel.qLHGR = 200e2;                                                      %INPUT

% initial fuel node radial thickness (m)
  fuel.dr0 = (fuel.rOut-fuel.rIn)/(fuel.nr-1);

% vector of fuel node radii (m)
  fuel.r0 = (fuel.rIn : fuel.dr0 : fuel.rOut)'; 

% vector of fuel node boundaries (m)
  fuel.r0_ = [fuel.rIn; interp1(fuel.r0,1.5:fuel.nr-0.5)'; fuel.rOut];

% vector of fuel node volumes (m3)
  fuel.v0 = pi * diff(fuel.r0_.^2) * g.h0; 

% vector of fuel node mass (kg)
  fuel.mass = fuel.rho * (1 - fuel.por) .* fuel.v0;

% power density (W/m3) assumed constant
  fuel.qV = fuel.qLHGR .* g.h0 ./ sum(fuel.v0);

% fission rate per unit volume of fuel (1/m3) assumed constant
  fuel.dFdt = fuel.qV / 1.60217656535e-13 / 200;

%-------------------------------------------------------------------------- 
% Initialize clad

% inner clad radius (m) as manufactured
  clad.rIn = 4.22e-3;                                                      %INPUT

% outer clad radius (m) as manufactured
  clad.rOut = 4.75e-3;                                                     %INPUT

% number of nodes in clad
  clad.nr = 5;                                                             %INPUT

% fast flux in cladding (1/cm2-s) assumed constant
  clad.fastFlux = 1e13;                                                    %INPUT

% outer cladding temperature (K) assumed constant
  clad.Tout = 600;                                                         %INPUT

% initial clad node radial thickness (m)
  clad.dr0 = (clad.rOut-clad.rIn)/(clad.nr-1);
  
% vector of clad node radii (m)
  clad.r0 = (clad.rIn : clad.dr0 : clad.rOut)';

% vector of clad node boundaries (m)
  clad.r0_ = [clad.rIn; interp1(clad.r0,1.5:clad.nr-0.5)'; clad.rOut];

% vector of clad node volumes (m3)
  clad.v0 = pi * diff(clad.r0_.^2) * g.h0;
  
%-------------------------------------------------------------------------- 
% Initialize gap

% fuel-clad effective roughness
  gap.rough = 6e-6;                                                        %INPUT

% gap width (m) as fabricated
  gap.dr0 = clad.rIn - fuel.rOut;

% gap status
  gap.open = 1;
  gap.clsd = 0;
  gap.depsh = 0;
  gap.depsz = 0;

% gap mid-radius (m) as fabricated
  gap.r0_ = (clad.rIn + fuel.rOut)/2;

%-------------------------------------------------------------------------- 
% Initialize inner gas

% fuel rod fission gas plenum temperature (K) assumed constant
  ingas.Tplenum = 600;                                                     %INPUT

% as-fabricated helium pressure inside fuel rod (MPa)
  p0 = 1;                                                                  %INPUT

% gas plenum volume (m3) assumed constant
  g.vGasPlenum = 10e-6;                                                    %INPUT

% gas gap volume (m3) as fabricated
  vGasGap = g.h0 * pi*(clad.rIn^2 - fuel.rOut.^2); 

% gas central void volume (m3) as fabricated
  vGasCentralVoid = g.h0 * pi*fuel.rIn.^2; 

% amount of helium inside fuel rod (mole) as-fabricated
  ingas.muHe0 = (p0*1e6) * (g.vGasPlenum + vGasGap + vGasCentralVoid) / (8.31 * 293);

%--------------------------------------------------------------------------
% ANONIMOUS FUNCTIONS
%--------------------------------------------------------------------------
  g.interpolateBetweenFuelNodes = @(y) interp1((1:fuel.nr)', y, (1.5:fuel.nr-0.5)');
  g.interpolateBetweenCladNodes = @(y) interp1((1:clad.nr)', y, (1.5:clad.nr-0.5)');
% Von Mises stress
  g.sigVM = @(sig) (0.5*(sig{1}-sig{2}).^2 + 0.5*(sig{1}-sig{3}).^2 + 0.5*(sig{2}-sig{3}).^2).^0.5; 

end