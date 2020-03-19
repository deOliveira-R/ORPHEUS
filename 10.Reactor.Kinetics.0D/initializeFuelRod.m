% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function initializes fuel rod parameters

function initializeFuelRod

% Global structures
  global g fuel clad gap ingas

%--------------------------------------------------------------------------
% FUEL ROD MATERIAL PROPERTIES
%--------------------------------------------------------------------------
% correlations for material properties (in folder 00.Lib)
  [fuel, gap, clad] = matpro(fuel, gap, clad);

%-------------------------------------------------------------------------- 
% Initialize fuel

% inner fuel radius (m) as fabricated
  fuel.rIn = 0;                                                            %INPUT

% outer fuel radius (m) as fabricated
  fuel.rOut = 4.12e-3;                                                     %INPUT

% number of radial nodes in fuel
  fuel.nr = 20;                                                            %INPUT

% initial porosity (-) assumed constant
  fuel.por = 0.05;                                                         %INPUT
  
% fuel burnup (MWd/kgHM)
  fuel.Bu = 0;                                                             %INPUT

% fission gas release (-) assumed constant
  fuel.FGR = 0.06;                                                         %INPUT

% initial fuel rod power: 2775 MW / 157 SAs / 256 fuel rods = 69140.625 W  %INPUT
  fuel.pow0 = 69140.625;

% fuel rod average linear heat generation rate (W/m)
  fuel.LHGR = fuel.pow0 / (g.dz0*g.nz);

% fuel rod axial power peaking = local LHGR / fuel.LHGR
  fuel.kz = ones(g.nz,1);                                                  %INPUT

% initial fuel node radial thickness (m)
  fuel.dr0 = (fuel.rOut-fuel.rIn)/(fuel.nr-1);

% vector of fuel node radii (m)
  fuel.r0 = (fuel.rIn : fuel.dr0 : fuel.rOut)'; 

% vector of fuel node boundaries (m)
  fuel.r0_ = [fuel.rIn; interp1(fuel.r0,1.5:fuel.nr-0.5)'; fuel.rOut];

% vector of fuel node volumes (m3)
  fuel.v0 = repmat(pi * diff(fuel.r0_.^2) * g.dz0,1,g.nz)';

% vector of fuel node mass (kg)
  fuel.mass = fuel.rho * (1 - fuel.por) .* fuel.v0;

% power density (W/m3)
  fuel.qV = repmat(fuel.LHGR * g.dz0 * fuel.kz ./ sum(fuel.v0,2),1,fuel.nr);

% fission rate per unit volume of fuel (1/m3) assumed constant
  fuel.dFdt = fuel.qV / 1.60217656535e-13 / 200;

% XS area of node boundary in radial direction (m2)
  fuel.a_ = repmat(2*pi * interp1(fuel.r0,1.5:fuel.nr-0.5)' .*  g.dz0, 1, g.nz);

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

% initial clad node radial thickness (m)
  clad.dr0 = (clad.rOut-clad.rIn)/(clad.nr-1);
  
% vector of clad node radii (m)
  clad.r0 = (clad.rIn : clad.dr0 : clad.rOut)';

% vector of clad node boundaries (m)
  clad.r0_ = [clad.rIn; interp1(clad.r0,1.5:clad.nr-0.5)'; clad.rOut];

% vector of clad node volumes (m3)
  clad.v0 = repmat(pi * diff(clad.r0_.^2) * g.dz0,1,g.nz)';

% vector of clad node mass (kg)
  clad.mass = clad.rho .* clad.v0;
  
% XS area of node boundary in radial direction (m2)
  clad.a_ = repmat(2*pi * interp1(clad.r0,1.5:clad.nr-0.5)' .*  g.dz0, 1, g.nz);

% clad failure logical flag
  clad.fail = false;

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

% XS area of node boundary in radial direction (m2)
  gap.a_ = repmat(2*pi * gap.r0_ .*  g.dz0, 1, g.nz);

%-------------------------------------------------------------------------- 
% Initialize inner gas

% as-fabricated helium pressure inside fuel rod (MPa)
  p0 = 1.2;                                                                %INPUT

% gas plenum volume (m3) assumed constant
  g.vGasPlenum = 10e-6;                                                    %INPUT

% gas gap volume (m3) as fabricated
  vGasGap = g.dz0 * g.nz * pi*(clad.rIn^2 - fuel.rOut.^2); 

% gas central void volume (m3) as fabricated
  vGasCentralVoid = g.dz0 * g.nz * pi*fuel.rIn.^2; 

% amount of helium inside fuel rod (mole) as-fabricated
  ingas.muHe0 = (p0*1e6) * (g.vGasPlenum + vGasGap + vGasCentralVoid) / (8.31 * 293);

%--------------------------------------------------------------------------
% ANONIMOUS FUNCTIONS
%--------------------------------------------------------------------------
  g.interpolateBetweenFuelNodes = @(y) interp2(1:fuel.nr,(1:g.nz)', y, 1.5:fuel.nr-0.5,(1:g.nz)');
  g.interpolateBetweenCladNodes = @(y) interp2(1:clad.nr,(1:g.nz)', y, 1.5:clad.nr-0.5,(1:g.nz)');
% Von Mises stress
  g.sigVM = @(sig) (0.5*(sig{1}-sig{2}).^2 + 0.5*(sig{1}-sig{3}).^2 + 0.5*(sig{2}-sig{3}).^2).^0.5; 

end