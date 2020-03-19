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
% reshape y to 2D temporal matrix
  y2D = reshape(y, g.nz, []);

% read fuel temperature (K)
  fuel.T = y2D(:, 1:fuel.nr); 
% cut off the temporal matrix
  y2D(:, 1:fuel.nr) = [];

% read clad temperature (K)
  clad.T = y2D(:, 1:clad.nr);
% cut off the temporal matrix
  y2D(:, 1:clad.nr) = [];

% read coolant enthalpy
  cool.h = y2D(:,1); 
% cut off the temporal matrix
  y2D(:,1) = [];

% read coolant pressure
  cool.p = y2D(:,1);
% cut off the temporal matrix
  y2D(:,1) = [];

% read clad plastic effective strain (-)
  clad.epsPeff = y2D(:,1:clad.nr); 
% cut off the temporal matrix
  y2D(:,1:clad.nr) = [];

  for i=1:3
    % read clad plastic strain components (-)
      clad.epsP{i}(:,1:clad.nr) = y2D(:,1:clad.nr); 
    % cut off the temporal vector
      y2D(:,1:clad.nr) = [];
  end

  for i=1:3
    % read clad stress components (MPa)
      clad.sig{i}(:,1:clad.nr) = y2D(:,1:clad.nr); 
    % cut off the temporal vector
      y2D(:,1:clad.nr) = [];
	% calculate clad stress components at the node boundaries
      clad.sig_{i} = g.interpolateBetweenCladNodes(clad.sig{i});
  end

%--------------------------------------------------------------------------
% THERMAL HYDRAULICS
%--------------------------------------------------------------------------
% mixture properties in nodes
  for i=1:g.nz
      [cool.pro(i),cool.Lsat(i),cool.Vsat(i)] = H2Oeq(cool.p(i), cool.h(i));
  end

% user-defined mixture velocity at inlet (m/s)
  cool.vel0 = interp1(cool.vel0t(1,:), cool.vel0t(2,:), t);

% user-defined coolant temperature at inlet (K)
  cool.T0 = interp1(cool.T0t(1,:), cool.T0t(2,:), t); 

% user-defined mixture pressure at outlet (MPs)
  cool.p0 = interp1(cool.p0t(1,:), cool.p0t(2,:), t); 

% mixture flowrate at inlet (kg/s) 
  cool.mdot0 = cool.vel0 * cool.aFlow0 * XSteam('rho_pT', cool.p(1)*10, cool.T0-273.15);

% mixture flowrate at junctions (kg/s) 
  cool.mdot_ = cool.mdot0 * ones(g.nz+1,1);

% mixture flowrate in nodes (kg/s) 
  cool.mdot = interp1(cool.mdot_, 1.5 : g.nz+0.5)';

% mixture velocity in nodes (m/s)
  cool.vel = cool.mdot ./ [cool.pro.rho]' ./ cool.aFlow;

% mixture reynolds number in nodes (-)
  cool.re = cool.dHyd .* abs(cool.vel) ./ [cool.pro.nu]';

%--------------------------------------------------------------------------
% PRESSURE CALCULATIONS
%--------------------------------------------------------------------------
% Churchill model for wall friction factor neglecting the wall roughness
  teta1 = -2.457*log( (7./cool.re).^0.9 );
  teta2 = 37530./cool.re;
  cool.fricFactorWall = 8*( (8./cool.re).^12 + 1 ./(teta1.^16 + teta2.^16).^1.5 ).^(1/12);
  
% wall friction pressure drop for mixture (MPa)
  cool.dpFricWall = cool.fricFactorWall .* g.dz0 ./ (2*cool.dHyd) .* [cool.pro.rho]' .* cool.vel.^2/2 * 1e-6;

% pressure drop in node due to gravity (MPa)
  cool.dpGrav = [cool.pro.rho]' * 9.81e-6 .* g.dz0;

% pressure drop in junctions due to acceleration (MPa)
  cool.dpAccel_ = -diff([cool.pro.rho]' .* cool.vel.^2) * 1e-6;

% pressure in junction (MPa)
  cool.p_ = cool.p0 + flip([0; cumsum(flip(cool.dpGrav + cool.dpFricWall))]);

% pressure in nodes (MPa)
  p_tilda = interp1(cool.p_, 1.5 : g.nz+0.5)' + flip([0; cumsum(flip(cool.dpAccel_))]);

% residual of pressure (Algebraic Equation)
  residualP = cool.p - p_tilda;

%--------------------------------------------------------------------------
% WALL HEAT EXCHANGE
%--------------------------------------------------------------------------
  Twall = clad.T(:,clad.nr);
  
% wall temperature minus mixture temperature (K)
  dTw = Twall - [cool.pro.T]';

% wall temperature minus saturation temperature (K)
  dTwSat = Twall - [cool.pro.Tsat]';

% difference in densities of liquid and gas at saturation (kg/m3)
  drhoSat = [cool.Lsat.rho]' - [cool.Vsat.rho]';

% latent heat of vapourization (J/kg)
  dhLat = [cool.Vsat.h]' - [cool.Lsat.h]';

% subcooling enthalpy (J/kg)
  dhSub = max([cool.Lsat.h]' - cool.h, 0);

% Prandtl number (-)
  cool.pr = [cool.pro.c_p]' .* [cool.pro.mu]' ./ [cool.pro.k]';

% subcooling correction multiplier (-)
  kSub = 1 + 0.1 * ([cool.Lsat.rho]' ./ [cool.Vsat.rho]').^0.75 .* dhSub ./ dhLat;

% void correction multiplier (-)
  kVoid = 1 - [cool.pro.void]';

% heat flux at single phase
  qw1PHASE =  max(4.36, 0.023 * cool.re.^0.8 .* cool.pr.^0.4) .* ([cool.pro.k]' ./ cool.dHyd) .* dTw;

% heat flux at nucleate boiling (Thom) (W/m2)
  qwNB = 2000 * abs(dTwSat) .* dTwSat .* exp(cool.p/4.34);

% critical heat flux (W/m2)
  qwCHF = 0.14 * dhLat .* (9.81 * [cool.Lsat.sig]' .* [cool.Vsat.rho]'.^2 .* drhoSat).^0.25 .* kSub .* kVoid;

% wall temperature at critical heat flux (K)
  cool.TCHF = [cool.pro.Tsat]' + sqrt(qwCHF .* exp(-cool.p/4.34) / 2000);

% film boiling heat flux (W/m2)
  qwFB = 0.25*(9.81 * [cool.Vsat.k]'.^2 .* [cool.Vsat.c_p]' .* drhoSat ./ [cool.Vsat.nu]').^0.333 .* kSub .* dTwSat;

%--------------------------------------------------------------------------
% FLOW REGIME SELECTION
%--------------------------------------------------------------------------
  for iz=1:g.nz
      if Twall(iz) < cool.Lsat(iz).T
         cool.regime(iz) = 1; % single-phase liquid forced convection
         cool.qw(iz) = qw1PHASE(iz);
         
      elseif Twall(iz) < cool.TCHF(iz)
         cool.regime(iz) = 2; % pre-CHF boiling
         cool.qw(iz) = max(qw1PHASE(iz),qwNB(iz));

      else
         cool.regime(iz) = 3; % post-CHF film boiling or single phase steam
         cool.qw(iz) = max(qw1PHASE(iz)*(1-kVoid(iz)),qwFB(iz)*kVoid(iz));
      end
  end

% gas generation rate in node (kg/s)
  cool.gamma = cool.qw' .* cool.areaHX ./ dhLat; 
  cool.gamma([cool.pro.x] >=1) = 0;
  cool.gamma([cool.pro.x] < 0) = 0;
  
% water enthalpy (J/kg) at core inlet
  cool.h0 = XSteam('h_pT', cool.p0*10, cool.T0-273.15)*1e3;

% specific enthalpy of mixture in junctions (upwind scheme) (J/kg)
  cool.h_ = [cool.h0; cool.h];

  rate_of_coolant_enthalpy = (-diff(cool.mdot_ .* cool.h_) + cool.qw' .* cool.areaHX) ./ ([cool.pro.rho]' .* cool.volFlow);

%--------------------------------------------------------------------------
% CALCULATION OF FUEL STRAIN (SIMPLIFIED, ONLY THERMAL EXPANSION)
%--------------------------------------------------------------------------
% fuel average temperature (K)
  fuel.Tavg = sum(fuel.v0 .* fuel.T,2)./sum(fuel.v0,2);

% fuel thermal expansion component (m/m)
  fuel.epsT = fuel.thExp(fuel.Tavg);

% fuel outer radius
  fuel.r = fuel.rOut*(1 + fuel.epsT);

% fuel node height
  fuel.dz = g.dz0*(1 + fuel.epsT);

%--------------------------------------------------------------------------
% CALCULATION OF CLAD STRAINS AND STRAIN RATES
%--------------------------------------------------------------------------
% clad average temperature (K)
  clad.Tavg = sum(clad.v0 .* clad.T,2)./sum(clad.v0,2);
  
% sum of stress components (MPa)
  sigSumC = clad.sig{1} + clad.sig{2} + clad.sig{3};

  for i=1:3
    % clad thermal expansion component (m/m)
      clad.epsT = clad.thExp(clad.T); 

	% clad elastic component (m/m)
      clad.epsE{i} = (clad.sig{i} - clad.nu*(sigSumC - clad.sig{i})) ./ clad.E(clad.T);

    % components of total clad deformation (m/m)
      clad.eps{i} = clad.epsT + clad.epsE{i} + clad.epsP{i};
      
      clad.eps_{i} = g.interpolateBetweenCladNodes(clad.eps{i});
  end
 
% Von Mises stress in clad
  clad.sigVM = g.sigVM(clad.sig) + 1e-6;

% effective plastic rate
  rate_of_eff_plastic_strain = ~clad.fail * min(0.1,1e-3*(clad.sigVM*1e6 ./ clad.K(clad.T) ./ abs(clad.epsPeff+1e-6).^clad.n(clad.T)).^(1./clad.m(clad.T)));

  for i=1:3
    % plastic components (m/m)
      rate_of_clad_plastic_strain{i} = rate_of_eff_plastic_strain .* (3/2) .* (clad.sig{i} - sigSumC/3) ./ clad.sigVM;
  end

%--------------------------------------------------------------------------
% UPDATE OF CLAD GEOMETRY
%--------------------------------------------------------------------------
% update node thicknesses (m)
  clad.dr = clad.dr0 * (1 + clad.eps_{1});

% update node radii (m)
  clad.r = repmat(clad.r0',g.nz,1) .* (1 + clad.eps{2});

% node inner boundary radii (m)
  clad.r_ = g.interpolateBetweenCladNodes(clad.r);

% node boundary radii (m)
  clad.r_plus = [clad.r(:,1), clad.r_, clad.r(:,clad.nr)];

% update node height
  clad.dz = repmat(g.dz0, g.nz, clad.nr) .* (1 + clad.eps{3});

% XS area of node boundary in radial direction (m2)
  clad.a_ = 2*pi * clad.r_ .*  g.interpolateBetweenCladNodes(clad.dz);

% node volume (m3)
  clad.v = pi*diff(clad.r_plus.^2,1,2) .* clad.dz;

%--------------------------------------------------------------------------
% UPDATE OF GAP GEOMETRY
%--------------------------------------------------------------------------
% gap width
  gap.dr = clad.r(:,1) - fuel.r;
   
% average gap radius (m)
  gap.r_ = (clad.r(1) + fuel.dz)/2;

% XS area of the mid-gap in radial direction (m2)
  gap.a_ = 2*pi * gap.r_ .* (clad.dz(1) + fuel.dz)/2;

% contact pressure
  gap.pContact = gap.clsd .* -clad.sig{1}(1);

%--------------------------------------------------------------------------
% CALCULATION OF INNER GAS VOLUME AND PRESSURE (SIMPLIFIED, NO CENTRAL
% VOLUME, NO FGR)
%--------------------------------------------------------------------------
% average gap temperature (K)
  gap.T = (fuel.T(fuel.nr) + clad.T(1))/2;

% gas gap volume (m3)
  vGasGap = fuel.dz * pi .* (clad.r(1).^2 - fuel.r.^2);

% gas plenum temperature (K)
  ingas.Tplenum = cool.T(g.nz);

% inner gas pressure (MPa)
  ingas.p = ingas.muHe0 * 8.31e-6 / ( g.vGasPlenum / ingas.Tplenum + ...
            + sum(vGasGap ./ gap.T) );

%--------------------------------------------------------------------------
% CLAD FAILURE
%--------------------------------------------------------------------------
% clad failure
  if clad.fail
     ingas.p = cool.p(end);
  end

% engineering stress in cladding
  clad.sigI = (ingas.p - cool.p(end)) * (clad.r(:,clad.nr) + clad.r(:,1))/2 ./ (clad.r(:,clad.nr) - clad.r(:,1));

%--------------------------------------------------------------------------
% RESIDUALS FOR ALGEBRAIC EQUATIONS FOR FUEL AND CLAD
%--------------------------------------------------------------------------
  residual_of_clad_stress_equilibrium_equation = diff(clad.sig{1},1,2) ./ clad.dr - (clad.sig_{2} - clad.sig_{1}) ./ clad.r_;
  residual_of_clad_strain_compatibility_equation = diff(clad.eps{2},1,2) ./ clad.dr - (clad.eps_{1} - clad.eps_{2}) ./ clad.r_;
  residual_of_ez_const_for_clad = diff(clad.eps{3},1,2);

%--------------------------------------------------------------------------
% RESIDUALS FOR ALGEBRAIC EQUATIONS FOR BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
  sigzCladIntegral = sum(clad.sig_{3} .* diff(clad.r.^2,1,2),2);
% radial stress equals minus coolant pressure on the clad outer surface
  residual_of_BC1 = clad.sig{1}(:,clad.nr) + cool.p;
% radial stress equals minus inner pressure on the clad inner surface
  residual_of_BC2 = clad.sig{1}(:,1) + ingas.p;
% balance of axial stress integral
  residual_of_BC3 = sigzCladIntegral - (ingas.p * clad.r(:,1).^2 - cool.p .* clad.r(:,clad.nr).^2);

%--------------------------------------------------------------------------
% FUEL ROD POWER CALCULATIONS
%--------------------------------------------------------------------------
% user-defined fuel rod power (W)
  fuel.pow = interp1(fuel.powt(1,:), fuel.powt(2,:), t); 

% fuel rod average linear heat generation rate (W/m)
  fuel.LHGR = fuel.pow / (g.dz0*g.nz);

% power density (W/m3)
  fuel.qV = repmat(fuel.LHGR * g.dz0 * fuel.kz ./ sum(fuel.v0,2),1,fuel.nr);

%--------------------------------------------------------------------------
% FUEL ROD THERMAL CALCULATIONS
%--------------------------------------------------------------------------
% temperature between fuel nodes (K)
  Tf_ = g.interpolateBetweenFuelNodes(fuel.T);

% temperature between clad nodes (K)
  Tc_ = g.interpolateBetweenCladNodes(clad.T);

% thermal conductivity of gas mixture (W/m-K), simplified, no fission gases
  gap.kGasMix = gap.kHe(gap.T);

% fuel-cladding gap conductance (W/m2-K) -- simplified, no dependence on contact pressure
  if gap.open
     gap.h = gap.kGasMix ./ gap.dr;
  else
     gap.h = gap.kGasMix ./ gap.rough;
  end

% radial heat flux at the mid-gap (W/m2)
  qGap_ = gap.h .* (fuel.T(:,fuel.nr) - clad.T(:,1));

% radial heat flux between fuel nodes (W/m2)
  qFuel_ = -fuel.k(Tf_, fuel.Bu, fuel.por) .* diff(fuel.T,1,2) ./ fuel.dr0; 

% radial heat transfer between fuel nodes with account for boundary conditions (W)
  QFuel_ = [zeros(g.nz,1), qFuel_.*fuel.a_', qGap_.*gap.a_];

% radial heat flux between clad nodes (W/m2)
  qClad_ = -clad.k(Tc_) .* diff(clad.T,1,2) ./ clad.dr0;

% radial heat transfer between clad nodes with account for boundary conditions (W)
  QClad_ = [qGap_.*gap.a_, qClad_.*clad.a_, cool.qw' .* cool.areaHX];

% time derivative of fuel temperature (K/s)
  rate_of_fuel_temperature = (-diff(QFuel_,1,2) + fuel.qV .* fuel.v0) ./ (fuel.mass .* fuel.cp(fuel.T));
  
% time derivative of clad temperature (K/s)
  rate_of_clad_temperature = -diff(QClad_,1,2) ./ (clad.rho .* clad.cp(clad.T) .* clad.v0);

  out = [...
         reshape(rate_of_fuel_temperature,[],1); ...
         reshape(rate_of_clad_temperature,[],1); ...
         rate_of_coolant_enthalpy; ...
         residualP; ...
         reshape(rate_of_eff_plastic_strain,[],1); ...
         reshape(rate_of_clad_plastic_strain{1},[],1); ...
         reshape(rate_of_clad_plastic_strain{2},[],1); ...
         reshape(rate_of_clad_plastic_strain{3},[],1); ...
		 reshape(residual_of_clad_stress_equilibrium_equation,[],1); ...
		 reshape(residual_of_clad_strain_compatibility_equation,[],1); ...
		 reshape(residual_of_ez_const_for_clad,[],1); ...
         residual_of_BC1; ...
         residual_of_BC2; ...
         residual_of_BC3; ...
        ];
end