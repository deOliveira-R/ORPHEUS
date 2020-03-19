% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates fuel behaviour under 1D approximation

function fuelBehaviour

  clear all;

% Global structures
  global fuel clad gap cool

%-------------------------------------------------------------------------- 
% Initialize fuel rod parameters
  initializeFuelRod;

%-------------------------------------------------------------------------- 
% Initialize coolant

% coolant pressure (MPa) assumed constant
  cool.p = 16;                                                             % INPUT

%--------------------------------------------------------------------------
% Time units:
  minute = 60; %s
  hour = 60*minute;
  day = 24*hour;
  year = 365*day;
  
%--------------------------------------------------------------------------
% Time vector:
  timeStep = 1*day;
  timeEnd = 6*year;

  timeVector = (0 :timeStep: timeEnd);

%--------------------------------------------------------------------------
% Diagonal "mass matrix" to multiply the time derivatives of the ODE
% system. Needed by ode15s to solve together ODEs (1) and AEs (0) 
  massMatrix = diag([...
                     ones(fuel.nr,1); ... % ODEs for fuel temperatures
                     ones(clad.nr-1,1); ... % ODEs for clad temperatures
                     1; ... % ODE for number of fissions per unit volume
                     ones(fuel.nr,1); ... % ODEs for volumetric fuel swelling
                     ones(fuel.nr*3,1); ... % ODEs for fuel creep strain components
                     ones(clad.nr*3,1); ... % ODEs for clad creep strain components
                     ones(clad.nr,1); ... % ODEs for clad plastic effective strain
                     ones(clad.nr*3,1); ... % ODEs for clad plastic strain components
                     zeros(fuel.nr*3+clad.nr*3,1); ... % AEs for fuel and clad stress components
                    ]);

%--------------------------------------------------------------------------
% ode15s options: function called after every successful integration step;
% mass matrix; events function; relative and absolute tolerances:
  options = odeset('OutputFcn',@(t,y,flag) writeResults(t,y,flag), ...
                   'Mass', massMatrix, ...
                   'Events',@gapClosureEvent, ...
                   'RelTol', 1e-6, ...                                     %INPUT
                   'AbsTol', 1e-4);                                        %INPUT

%--------------------------------------------------------------------------
% Initial guess for unknowns
  y = [... 
       clad.Tout *ones(fuel.nr,1); ... % fuel temperatures
       clad.Tout *ones(clad.nr-1,1); ... % clad temperatures
       0; ... % number of fissions per unit volume
       zeros(fuel.nr,1); ... % volumetric fuel swelling
       zeros(fuel.nr*3,1); ... % fuel creep strain components
       zeros(clad.nr*3,1); ... % clad creep strain components
       zeros(clad.nr,1); ... % clad plastic effective strain
       zeros(clad.nr*3,1); ... % clad plastic strain components
       zeros((fuel.nr+clad.nr)*3,1); ... % fuel and clad stress components
      ]; 

%--------------------------------------------------------------------------
% ode15s solves the system of ordinary differential and algebraic
% equations, using funRHS for calculating the right-hand side vector;
% timeVector for outputing the results, vector y as an initial guess; and
% options (see above)
  [timeEnd1,y] = ode15s(@(t,x) funRHS(t,x), timeVector, y, options);

  if timeEnd1 < timeEnd
   % this is the gap closure... go on with the closed gap regime
     fprintf('Fuel-clad gap closed...\n');
     gap.open = 0;
     gap.clsd = 1;
     gap.depsz = clad.eps{3}(1) - fuel.eps{3}(fuel.nr);
     gap.depsh = clad.eps{2}(1) - fuel.eps{2}(fuel.nr);

     options = odeset('OutputFcn',@(t,y,flag) writeResults(t,y,flag), ...
                      'Mass', massMatrix, ...
                      'RelTol', 1e-4, ...                                  %INPUT
                      'AbsTol', 1e-4);                                     %INPUT

     ode15s(@(t,x) funRHS(t,x), (timeEnd1(end,1) :timeStep: timeEnd), y(end,:)', options);
  end

  fprintf('Plot the results...\n');

%--------------------------------------------------------------------------
% Plot the results
  s = results;

  f = figure('visible','off');
  plot(s.fuel.r(:,1),s.fuel.T(:,1),'-ob', ...
       s.fuel.r(:,end),s.fuel.T(:,end),'-or', ...
       s.clad.r(:,1),s.clad.T(:,1),'-ob', ...
       s.clad.r(:,end),s.clad.T(:,end),'-or');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Temperature (C)');
  legend('at BOC','at EOC','Location','best');
  saveas(f, 'Fig_radial_profiles_of_temperatures.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,1),s.fuel.sigr(:,1),'-or', ...
       s.fuel.r(:,1),s.fuel.sigh(:,1),'-og', ...
       s.fuel.r(:,1),s.fuel.sigz(:,1),'-ob', ...
       s.fuel.r(:,1),s.fuel.sigVM(:,1),'-ok', ...
       s.clad.r(:,1),s.clad.sigVM(:,1),'-oK', ...
       s.clad.r(:,1),s.clad.sigr(:,1),'-or', ...
       s.clad.r(:,1),s.clad.sigh(:,1),'-og', ...
       s.clad.r(:,1),s.clad.sigz(:,1),'-ob')
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Stress (MPa)');
  legend('radial','hoop','axial','Von Misses', ...
         'Location','best');
  saveas(f, 'Fig_radial_profiles_of_stresses@BOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,end),s.fuel.sigr(:,end),'-or', ...
       s.fuel.r(:,end),s.fuel.sigh(:,end),'-og', ...
       s.fuel.r(:,end),s.fuel.sigz(:,end),'-ob', ...
       s.fuel.r(:,end),s.fuel.sigVM(:,end),'-ok', ...
       s.clad.r(:,end),s.clad.sigVM(:,end),'-oK', ...
       s.clad.r(:,end),s.clad.sigr(:,end),'-or', ...
       s.clad.r(:,end),s.clad.sigh(:,end),'-og', ...
       s.clad.r(:,end),s.clad.sigz(:,end),'-ob')
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Stress (MPa)');
  legend('radial','hoop','axial','Von Misses', ...
         'Location','best');
  saveas(f, 'Fig_radial_profiles_of_stresses@EOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,1),s.fuel.epsT(:,1),'-.k', ...
       s.fuel.r(:,1),s.fuel.epsS(:,1)/3,'-k', ...
       s.fuel.r(:,1),s.fuel.epsrE(:,1),'-.r', ...
       s.fuel.r(:,1),s.fuel.epsrC(:,1),'-r', ...
       s.fuel.r(:,1),s.fuel.epsr(:,1),'-or', ...
       s.clad.r(:,1),s.clad.epsT(:,1),'-.k', ...
       s.clad.r(:,1),s.clad.epsrE(:,1),'-.r', ...
       s.clad.r(:,1),s.clad.epsrC(:,1),'-r', ...
       s.clad.r(:,1),s.clad.epsr(:,1),'-or');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Radial strain (%)');
  legend('linear thermal', 'linear swelling', 'elastic', 'creep', 'total', 'Location','best');
  saveas(f, 'Fig_radial_profiles_of_radial_strains@BOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,end),s.fuel.epsT(:,end),'-.k', ...
       s.fuel.r(:,end),s.fuel.epsS(:,end)/3,'-k', ...
       s.fuel.r(:,end),s.fuel.epsrE(:,end),'-.r', ...
       s.fuel.r(:,end),s.fuel.epsrC(:,end),'-r', ...
       s.fuel.r(:,end),s.fuel.epsr(:,end),'-or', ...
       s.clad.r(:,end),s.clad.epsT(:,end),'-.k', ...
       s.clad.r(:,end),s.clad.epsrE(:,end),'-.r', ...
       s.clad.r(:,end),s.clad.epsrC(:,end),'-r', ...
       s.clad.r(:,end),s.clad.epsr(:,end),'-or');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Radial strain (%)');
  legend('linear thermal', 'linear swelling', 'elastic', 'creep', 'total', 'Location','best');
  saveas(f, 'Fig_radial_profiles_of_radial_strains@EOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,1),s.fuel.epsT(:,1),'-.k', ...
       s.fuel.r(:,1),s.fuel.epsS(:,1)/3,'-k', ...
       s.fuel.r(:,1),s.fuel.epshE(:,1),'-.g', ...
       s.fuel.r(:,1),s.fuel.epshC(:,1),'-g', ...
       s.fuel.r(:,1),s.fuel.epsh(:,1),'-og', ...
       s.clad.r(:,1),s.clad.epsT(:,1),'-.k', ...
       s.clad.r(:,1),s.clad.epshE(:,1),'-.g', ...
       s.clad.r(:,1),s.clad.epshC(:,1),'-g', ...
       s.clad.r(:,1),s.clad.epsh(:,1),'-og');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Hoop strain (%)');
  legend('linear thermal', 'linear swelling', 'elastic', 'creep', 'total', 'Location','best');
  saveas(f, 'Fig_radial_profiles_of_hoop_strains@BOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,end),s.fuel.epsT(:,end),'-.k', ...
       s.fuel.r(:,end),s.fuel.epsS(:,end)/3,'-k', ...
       s.fuel.r(:,end),s.fuel.epshE(:,end),'-.g', ...
       s.fuel.r(:,end),s.fuel.epshC(:,end),'-g', ...
       s.fuel.r(:,end),s.fuel.epsh(:,end),'-og', ...
       s.clad.r(:,end),s.clad.epsT(:,end),'-.k', ...
       s.clad.r(:,end),s.clad.epshE(:,end),'-.g', ...
       s.clad.r(:,end),s.clad.epshC(:,end),'-g', ...
       s.clad.r(:,end),s.clad.epsh(:,end),'-og');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Hoop strain (%)');
  legend('linear thermal', 'linear swelling', 'elastic', 'creep', 'total', 'Location','best');
  saveas(f, 'Fig_radial_profiles_of_hoop_strains@EOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,1),s.fuel.epsT(:,1),'-.k', ...
       s.fuel.r(:,1),s.fuel.epsS(:,1)/3,'-k', ...
       s.fuel.r(:,1),s.fuel.epszE(:,1),'-.b', ...
       s.fuel.r(:,1),s.fuel.epszC(:,1),'-b', ...
       s.fuel.r(:,1),s.fuel.epsz(:,1),'-ob', ...
       s.clad.r(:,1),s.clad.epsT(:,1),'-.k', ...
       s.clad.r(:,1),s.clad.epszE(:,1),'-.b', ...
       s.clad.r(:,1),s.clad.epszC(:,1),'-b', ...
       s.clad.r(:,1),s.clad.epsz(:,1),'-ob');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Axial strain (%)');
  legend('linear thermal', 'linear swelling', 'elastic', 'creep', 'total', 'Location','best');
  saveas(f, 'Fig_radial_profiles_of_axial_strains@BOC.pdf');

  f = figure('visible','off');
  plot(s.fuel.r(:,end),s.fuel.epsT(:,end),'-.k', ...
       s.fuel.r(:,end),s.fuel.epsS(:,end)/3,'-k', ...
       s.fuel.r(:,end),s.fuel.epszE(:,end),'-.b', ...
       s.fuel.r(:,end),s.fuel.epszC(:,end),'-b', ...
       s.fuel.r(:,end),s.fuel.epsz(:,end),'-ob', ...
       s.clad.r(:,end),s.clad.epsT(:,end),'-.k', ...
       s.clad.r(:,end),s.clad.epszE(:,end),'-.b', ...
       s.clad.r(:,end),s.clad.epszC(:,end),'-b', ...
       s.clad.r(:,end),s.clad.epsz(:,end),'-ob');
  ylim(ylim);
  grid on;
  xlabel('Radius (mm)');
  ylabel('Axial strain (%)');
  legend('linear thermal', 'linear swelling', 'elastic', 'creep', 'total', 'Location','best');
  saveas(f, 'Fig_radial_profiles_of_axial_strains@EOC.pdf');

  f = figure('visible','off');
  plot(s.timeY,s.fuel.r(end,:),'-r', ...
       s.timeY,s.clad.r(1,:),'-b', ...
       s.timeY,s.clad.r(end,:),'--b')
  ylim(ylim);
  grid on;
  xlabel('Time (years)');
  ylabel('Radius (mm)');
  legend('fuel outer surface','clad inner surface','clad outer surface', ...
         'Location','best');
  saveas(f, 'Fig_fuel&clad_radii_vs_time.pdf');

  f = figure('visible','off');
  plot(s.timeY,s.ingas.p(end,:),'-r');
  ylim(ylim);
  grid on;
  xlabel('Time (years)');
  ylabel('Inner gas pressure (MPa)');
  saveas(f, 'Fig_inner_gas_pressure_vs_time.pdf');

  f = figure('visible','off');
  plot(s.timeY,s.fuel.T(1,:),'-r', ...
       s.timeY,s.fuel.T(end,:),'--r', ...
       s.timeY,s.clad.T(1,:),'-b',...
       s.timeY,s.clad.T(end,:),'--b')
  ylim(ylim);
  grid on;
  xlabel('Time (years)');
  ylabel('Temperature (C)');
  legend('fuel inner surface','fuel outer surface','clad inner surface', 'clad outer surface', ...
         'Location','best');
  saveas(f, 'Fig_fuel&clad_temperature_vs_time.pdf');

  f = figure('visible','off');
  plot(s.timeY,s.gap.dr(:),'-r')
  ylim(ylim);
  grid on;
  xlabel('Time (years)');
  ylabel('Gap width (mkm)');
  saveas(f, 'Fig_gap_width_vs_time.pdf');

  f = figure('visible','off');
  plot(s.timeY,s.gap.pcontact(:),'-b')
  ylim(ylim);
  grid on;
  xlabel('Time (years)');
  ylabel('Contact pressure (MPa)');
  saveas(f, 'Fig_contact_pressure_vs_time.pdf');

  f = figure('visible','off');
  plot(s.timeY,s.fuel.dz(:),'-r', ...
       s.timeY,s.clad.dz(:),'-b');
  ylim(ylim);
  grid on;
  xlabel('Time (years)');
  ylabel('Height (mm)');
  legend('fuel','clad', ...
         'Location','best');
  saveas(f, 'Fig_fuel&clad_height_vs_time.pdf');

  fprintf('Done.\n');

end