% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function is called by ode15s and, given the time t and the vector of
% uknowns y, returns the right-hand sides of the matrix

function out = funRHS(t, y)

  global g fuel gap ingas clad cool

%--------------------------------------------------------------------------
% READ VECTOR OF UNKNOWNS
%--------------------------------------------------------------------------
% copy y to temporal vector
  y_ = y;

% read fuel temperature (K)
  fuel.T = y_(1:fuel.nr); 
% cut off the temporal vector
  y_(1:fuel.nr) = [];

% read clad temperature (K)
  clad.T = [y_(1:clad.nr-1); clad.Tout];  
% cut off the temporal vector
  y_(1:clad.nr-1) = [];
  
% read number of fissions per unit volume (1/m3)
  fuel.F = y_(1)*1e27;
% cut off the temporal vector
  y_(1) = [];

% read fuel swelling (-)
  fuel.epsS = y_(1:fuel.nr)/100;
% cut off the temporal vector
  y_(1:fuel.nr) = [];  

  for i=1:3
    % read fuel creep strain components (-)
      fuel.epsC{i} = y_(1:fuel.nr)/100;
    % cut off the temporal vector
      y_(1:fuel.nr) = [];
  end

  for i=1:3
    % read clad creep strain components (-)
      clad.epsC{i} = y_(1:clad.nr)/100; 
    % cut off the temporal vector
      y_(1:clad.nr) = [];
  end

% read clad plastic effective strain (-)
  clad.epsPeff = y_(1:clad.nr)/100; 
% cut off the temporal vector
  y_(1:clad.nr) = [];

  for i=1:3
    % read clad plastic strain components (-)
      clad.epsP{i} = y_(1:clad.nr)/100; 
    % cut off the temporal vector
      y_(1:clad.nr) = [];
  end

  for i=1:3
    % read fuel stress components (MPa)
      fuel.sig{i} = y_(1:fuel.nr); 
    % cut off the temporal vector
      y_(1:fuel.nr) = [];
	% calculate fuel stress components at the node boundaries
      fuel.sig_{i} = g.interpolateBetweenFuelNodes(fuel.sig{i});
  end

  for i=1:3
    % read clad stress components (MPa)
      clad.sig{i} = y_(1:clad.nr); 
    % cut off the temporal vector
      y_(1:clad.nr) = [];  
	% calculate clad stress components at the node boundaries
      clad.sig_{i} = g.interpolateBetweenCladNodes(clad.sig{i});
  end

%--------------------------------------------------------------------------
% FISSION RATE (1/s)AND FUEL BURNUP (MWd/kgUO2)
%--------------------------------------------------------------------------
  rate_of_fission = fuel.dFdt;
  fuel.Bu = fuel.F * 200 * 1.60217656535e-13 / 8.64e10 / fuel.rho; 
  
%--------------------------------------------------------------------------
% CALCULATION OF FUEL STRAINS AND STRAIN RATES
%--------------------------------------------------------------------------
% sum of fuel stress components (MPa)
  sigSumF = fuel.sig{1} + fuel.sig{2} + fuel.sig{3}; 

  for i=1:3
    % fuel thermal expansion component (m/m)
      fuel.epsT = fuel.thExp(fuel.T);

    % fuel elastic components (m/m)
      fuel.epsE{i} = (fuel.sig{i} - fuel.nu*(sigSumF - fuel.sig{i})) ./ fuel.E(fuel.T,fuel.por);

    % components of total fuel deformation (m/m)
      fuel.eps{i} = fuel.epsT + fuel.epsE{i} + fuel.epsC{i} + fuel.epsS/3; 

      fuel.eps_{i} = g.interpolateBetweenFuelNodes(fuel.eps{i});
  end

% fuel swelling rate (1/s)
  rate_of_fuel_swelling = fuel.swelRate(rate_of_fission, fuel.F, fuel.T);

% Von Mises stress in fuel
  fuel.sigVM = g.sigVM(fuel.sig) + 1e-6; 

% effective fuel creep rate
  effCreepRate = fuel.creepRate(fuel.sigVM, fuel.T);
  for i=1:3
    % fuel creep rate components (1/s) according to Prandtl-Reuss flow rule
      rate_of_fuel_creep{i} = effCreepRate .* (fuel.sig{i} - sigSumF/3) ./ fuel.sigVM; 
  end

%--------------------------------------------------------------------------
% CALCULATION OF CLAD STRAINS AND STRAIN RATES
%--------------------------------------------------------------------------
% sum of stress components (MPa)
  sigSumC = clad.sig{1} + clad.sig{2} + clad.sig{3};

  for i=1:3
    % clad thermal expansion component (m/m)
      clad.epsT = clad.thExp(clad.T); 

	% clad elastic component (m/m)
      clad.epsE{i} = (clad.sig{i} - clad.nu*(sigSumC - clad.sig{i})) ./ clad.E(clad.T);

    % components of total clad deformation (m/m)
      clad.eps{i} = clad.epsT + clad.epsE{i} + clad.epsP{i} + clad.epsC{i};
      
      clad.eps_{i} = g.interpolateBetweenCladNodes(clad.eps{i});
  end

% Von Mises stress in clad
  clad.sigVM = g.sigVM(clad.sig) + 1e-6;

% effective creep rate
  effCreepRate = clad.creepRate(clad.sigVM, clad.T);

% effective plastic rate
  rate_of_eff_plastic_strain = 1e-3*(clad.sigVM*1e6 ./ clad.K(clad.T) ./ abs(clad.epsPeff+1e-6).^clad.n(clad.T)).^(1./clad.m(clad.T));

  for i=1:3
    % creep components (m/m)
      rate_of_clad_creep{i} = effCreepRate .* (3/2) .* (clad.sig{i} - sigSumC/3) ./ clad.sigVM;

    % plastic components (m/m)
      rate_of_clad_plastic_strain{i} = (rate_of_eff_plastic_strain >= clad.sigVM./clad.E(clad.T)) .* ...
                                        rate_of_eff_plastic_strain .* (3/2) .* (clad.sig{i} - sigSumC/3) ./ clad.sigVM;
  end

%--------------------------------------------------------------------------
% UPDATE OF FUEL AND CLAD GEOMETRY
%--------------------------------------------------------------------------
% update node thicknesses (m)
  fuel.dr = fuel.dr0 * (1 + fuel.eps_{1});
  clad.dr = clad.dr0 * (1 + clad.eps_{1});

% update node radii (m)
  fuel.r = fuel.r0 .* (1 + fuel.eps{2});
  clad.r = clad.r0 .* (1 + clad.eps{2});

% node inner boundary radii (m)
  fuel.r_ = g.interpolateBetweenFuelNodes(fuel.r);
  clad.r_ = g.interpolateBetweenCladNodes(clad.r);

% node boundary radii (m)
  fuel.r_plus = [fuel.r(1); fuel.r_; fuel.r(fuel.nr)];
  clad.r_plus = [clad.r(1); clad.r_; clad.r(clad.nr)];

% update node height
  fuel.dz = g.h0 * (1 + fuel.eps{3});
  clad.dz = g.h0 * (1 + clad.eps{3});

% XS area of node boundary in radial direction (m2)
  fuel.a_ = 2*pi * fuel.r_ .*  g.interpolateBetweenFuelNodes(fuel.dz);
  clad.a_ = 2*pi * clad.r_ .*  g.interpolateBetweenCladNodes(clad.dz);

% node volume (m3)
  fuel.v = pi*diff(fuel.r_plus.^2) .* fuel.dz;
  clad.v = pi*diff(clad.r_plus.^2) .* clad.dz;

%--------------------------------------------------------------------------
% UPDATE OF GAP GEOMETRY
%--------------------------------------------------------------------------
% gap width
  gap.dr = clad.r(1) - fuel.r(fuel.nr);
   
% average gap radius (m)
  gap.r_ = (clad.r(1) + fuel.r(fuel.nr))/2;

% XS area of the mid-gap in radial direction (m2)
  gap.a_ = 2*pi * gap.r_ .* (clad.dz(1) + fuel.dz(fuel.nr))/2;

% contact pressure
  gap.pContact = gap.clsd .* -clad.sig{1}(1);

%--------------------------------------------------------------------------
% CALCULATION OF INNER GAS VOLUME, COMPOSITION AND PRESSURE
%--------------------------------------------------------------------------
% average gap temperature (K)
  gap.T = (fuel.T(fuel.nr) + clad.T(1))/2;

% gap temperature drop (K)
  dTGasGap = fuel.T(fuel.nr) - clad.T(1);

% gas gap volume (m3)
  vGasGap = fuel.dz(fuel.nr) * pi .* (clad.r(1).^2 - fuel.r(fuel.nr).^2);

% gas central void volume (m3)
  vGasCentralVoid = fuel.dz(1) * pi .* fuel.r(1).^2;

% fraction of fission gases He, Kr, Xe in total number of fission products
  fracFG = [0.01; 0.02; 0.23];

% [He; Kr; Xe] generated (mol)
  fuel.muGen = fuel.F * sum(fuel.v0) / 6.022e23 * fracFG;

% [He; Kr; Xe] released (mol)
  fuel.muRel = fuel.muGen * fuel.FGR;

% [He; Kr; Xe] amount (mol)
  ingas.mu = fuel.muRel + [1; 0; 0] .* ingas.muHe0;

% inner gas pressure (MPa)
  ingas.p = sum(ingas.mu) * 8.31e-6 / ( g.vGasPlenum / ingas.Tplenum + ...
            sum(vGasCentralVoid ./ fuel.T(1) + vGasGap ./ gap.T) );

%--------------------------------------------------------------------------
% RESIDUALS FOR ALGEBRAIC EQUATIONS FOR FUEL AND CLAD
%--------------------------------------------------------------------------
  residual_of_fuel_stress_equilibrium_equation = diff(fuel.sig{1}) ./ fuel.dr - (fuel.sig_{2} - fuel.sig_{1}) ./ fuel.r_;
  residual_of_fuel_strain_compatibility_equation = diff(fuel.eps{2}) ./ fuel.dr - (fuel.eps_{1} - fuel.eps_{2}) ./ fuel.r_;
  residual_of_ez_const_for_fuel = diff(fuel.eps{3});

  residual_of_clad_stress_equilibrium_equation = diff(clad.sig{1}) ./ clad.dr - (clad.sig_{2} - clad.sig_{1}) ./ clad.r_;
  residual_of_clad_strain_compatibility_equation = diff(clad.eps{2}) ./ clad.dr - (clad.eps_{1} - clad.eps_{2}) ./ clad.r_;
  residual_of_ez_const_for_clad = diff(clad.eps{3});

%--------------------------------------------------------------------------
% RESIDUALS FOR ALGEBRAIC EQUATIONS FOR BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
  sigzFuelIntegral = sum(fuel.sig_{3} .* diff(fuel.r.^2));
  sigzCladIntegral = sum(clad.sig_{3} .* diff(clad.r.^2));
  
% radial stress equals minus inner pressure (there is a central void) or
% radial and hoop stress components are equal (no central void)
  residual_of_BC1 = (fuel.r(1) >   0) .* (fuel.sig{1}(1) + ingas.p) + ... 
                    (fuel.r(1) ==  0) .* (fuel.sig{1}(1) - fuel.sig{2}(1));

% radial stress equals minus coolant pressure on the outer surface
  residual_of_BC2 = clad.sig{1}(clad.nr) + cool.p;

% radial stress at outer fuel surface equals minus inner pressure (gap is
% open) or clad inner radius equals fuel outer radius (gap is closed)
  residual_of_BC3 = gap.open * (fuel.sig{1}(fuel.nr) + ingas.p) + ...
                    gap.clsd * ((clad.sig{1}(1) - fuel.sig{1}(fuel.nr)) ./ gap.dr - ((fuel.sig{2}(fuel.nr) + clad.sig{2}(1))/2 - (fuel.sig{1}(fuel.nr) + clad.sig{1}(1))/2) ./ gap.r_);

% radial stress at inner clad surface equals minus inner pressure (gap is
% open) or clad inner radius equals fuel outer radius (gap is closed)
  residual_of_BC4 = gap.open * (clad.sig{1}(1) + ingas.p) + ...
                    gap.clsd * ((clad.eps{2}(1) - fuel.eps{2}(fuel.nr) - gap.depsh) ./ gap.dr - ((fuel.eps{1}(fuel.nr) + clad.eps{1}(1))/2 - (fuel.eps{2}(fuel.nr) + clad.eps{2}(1))/2) ./ gap.r_);

% integral of axial stress in fuel = 0 (open gap) or increment of axial strain are
% equal in contact (closed gap)
  residual_of_BC5 = gap.open * (sigzFuelIntegral) + ...
                    gap.clsd * (clad.eps{3}(1) - fuel.eps{3}(fuel.nr) - gap.depsz);

% integral of axial stress in cladding = pressure differential
% (open gap) or hoop stresses are equal in contact (closed gap)
  residual_of_BC6 = gap.open * (sigzCladIntegral - (ingas.p * clad.r(1).^2 - cool.p .* clad.r(clad.nr).^2)) + ...
                    gap.clsd * (sigzFuelIntegral + sigzCladIntegral - (ingas.p * clad.r(1).^2 - cool.p .* clad.r(clad.nr).^2));
                
%--------------------------------------------------------------------------
% THERMAL CALCULATIONS
%--------------------------------------------------------------------------
% temperature between fuel nodes (K)
  Tf_ = g.interpolateBetweenFuelNodes(fuel.T);

% temperature between clad nodes (K)
  Tc_ = g.interpolateBetweenCladNodes(clad.T);

% thermal conductivity of gas mixture (W/m-K)
  gap.kGasMix = gap.kGasMixFun({gap.kHe(gap.T); gap.kXe(gap.T); gap.kKr(gap.T)}, ...
                               [ingas.mu(1);    ingas.mu(2);    ingas.mu(3)]/sum(ingas.mu), ...
                               [2;              36;             54]);

% fuel-cladding gap conductance (W/m2-K) -- simplified, no dependence on contact pressure
  if gap.open
     gap.h = gap.kGasMix ./ gap.dr;
  else
     gap.h = gap.kGasMix ./ gap.rough;
  end

% radial heat flux at the mid-gap (W/m2)
  qGap_ = gap.h .* dTGasGap;

% radial heat flux between fuel nodes (W/m2)
  qFuel_ = -fuel.k(Tf_, fuel.Bu, fuel.por) .* diff(fuel.T) ./ fuel.dr; 

% radial heat transfer between fuel nodes with account for boundary conditions (W)  
  QFuel_ = [0; qFuel_.*fuel.a_; qGap_.*gap.a_];

% radial heat flux between clad nodes (W/m2)
  qClad_ = -clad.k(Tc_) .* diff(clad.T) ./ clad.dr;

% radial heat transfer between clad nodes with account for boundary conditions (W)
  QClad_ = [qGap_.*gap.a_; qClad_.*clad.a_; 0];

% time derivative of fuel temperature (K/s)
  rate_of_fuel_temperature = (-diff(QFuel_) + fuel.qV .* fuel.v0) ./ (fuel.mass .* fuel.cp(fuel.T));
  
% time derivative of clad temperature (K/s)
  rate_of_clad_temperature = -diff(QClad_) ./ (clad.rho .* clad.cp(clad.T) .* clad.v0);
  
% outer surface clad temperature specified as boundary condition: exclude it
  rate_of_clad_temperature(clad.nr) = [];

%--------------------------------------------------------------------------
% CONSTRUCTION OF OUTPUT VECTOR (1D)
%--------------------------------------------------------------------------
  out = [ ...
          rate_of_fuel_temperature; ...
          rate_of_clad_temperature; ...
          rate_of_fission/1e27; ...
          rate_of_fuel_swelling*100; ...
          rate_of_fuel_creep{1}*100; ...
          rate_of_fuel_creep{2}*100; ...
          rate_of_fuel_creep{3}*100; ...
          rate_of_clad_creep{1}*100; ...
          rate_of_clad_creep{2}*100; ...
          rate_of_clad_creep{3}*100; ...
		  rate_of_eff_plastic_strain*100; ...
          rate_of_clad_plastic_strain{1}*100; ...
          rate_of_clad_plastic_strain{2}*100; ...
          rate_of_clad_plastic_strain{3}*100; ...

          residual_of_fuel_stress_equilibrium_equation; ...
          residual_of_fuel_strain_compatibility_equation; ...
          residual_of_ez_const_for_fuel; ...
          residual_of_clad_stress_equilibrium_equation; ...
          residual_of_clad_strain_compatibility_equation; ...
          residual_of_ez_const_for_clad; ...          
          residual_of_BC1; ...
          residual_of_BC2; ...
          residual_of_BC3; ...
          residual_of_BC4; ...
          residual_of_BC5; ...
          residual_of_BC6; ...
        ];
end