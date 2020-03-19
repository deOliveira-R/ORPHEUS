% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates coupled 0D reactor kinetics, 1D thermal hydraulics 
% and fuel behaviour for one PWR-like subassembly under RIA conditions

function reactorKinetics

  clear all;

% Global structures
  global g fuel gap clad cool kin

% Start stopwatch timer
  g.timerValue = tic;

%-------------------------------------------------------------------------- 
% Initialize axial mesh

% number of axial nodes
  g.nz = 2;

% initial height of the axial node (m)
  g.dz0 = 1.5;

%-------------------------------------------------------------------------- 
% Initialize fuel rod parameters
  initializeFuelRod;

%-------------------------------------------------------------------------- 
% Initialize coolant parameters
  initializeCoolant;

%-------------------------------------------------------------------------- 
% Initialize neutronics parameters
  initializeNeutronics;

%--------------------------------------------------------------------------
% Time vector for steady state
  g.timeStep = 1; %s
  g.timeEnd1 = 100;  %s
  g.timeVector = (0 :g.timeStep: g.timeEnd1);
  
  g.tran = false;
  
%--------------------------------------------------------------------------
% Diagonal "mass matrix" to multiply the time derivatives of the ODE
% system. Needed by ode15s to solve together ODEs (1) and AEs (0) 
  massMatrix = diag([... 
                     1; ... % ODE for power
                     ones(6,1); ... % ODEs for delayed neutron precursors density
                     ones(g.nz*fuel.nr,1); ... % ODEs for fuel temperatures
                     ones(g.nz*clad.nr,1); ... % ODEs for clad temperatures
                     ones(g.nz,1); ... % ODEs for coolant enthalpy
                     zeros(g.nz,1); ...  % AEs for coolant pressures
                     ones(g.nz*clad.nr,1); ... % ODEs for clad plastic effective strain
                     ones(g.nz*clad.nr*3,1); ... % ODEs for clad plastic strain components
                     zeros(g.nz*clad.nr*3,1); ... % AEs for clad stress components
                    ]);

%--------------------------------------------------------------------------
% ode15s options: function called after every successful integration step;
% mass matrix; statistics flag; relative and absolute tolerances:
  options = odeset('OutputFcn',@(t,y,flag) writeResults(t,y,flag), ...
                   'Mass', massMatrix, ...
                   'MaxStep', 10, ...
                   'RelTol', 1e-6, ...                                     %INPUT
                   'AbsTol', 1e-4);                                        %INPUT
  
%--------------------------------------------------------------------------
% Initial guess for unknowns
  y = [...
       fuel.pow0; ...
       kin.cDNP'; ...
       cool.T0 *ones(g.nz*fuel.nr,1); ... % fuel temperatures
       cool.T0 *ones(g.nz*clad.nr,1); ... % fuel and clad temperatures
       cool.h; ... % coolant enthalpy
       cool.p; ... % coolant pressure
       zeros(g.nz*clad.nr,1); ... % clad plastic effective strain
       zeros(g.nz*clad.nr*3,1); ... % clad plastic strain components
       zeros(g.nz*clad.nr*3,1); ... % clad stress components
      ];

%--------------------------------------------------------------------------
% ode15s solves the system of ordinary differential and algebraic
% equations, using funRHS for calculating the right-hand side vector;
% timeVector for outputing the results, vector y as an initial guess; and
% options (see above)
  [~,y] = ode15s(@(t,x) funRHS(t,x), g.timeVector, y, options);

%--------------------------------------------------------------------------
% Time vector for transient before gap closure
  g.timeStep = 0.1; %s
  g.timeEnd2 = 120;  %s                                                    % INPUT
  g.timeVector = (g.timeEnd1 :g.timeStep: g.timeEnd2);
  
  g.tran = true;

  fprintf('transient starts\n');
  options = odeset('OutputFcn',@(t,y,flag) writeResults(t,y,flag), ...
                   'Mass', massMatrix, ...
                   'Events',@gapClosureEvent, ...
                   'MaxStep', 1e-3, ...
                   'RelTol', 1e-6, ...                                     %INPUT
                   'AbsTol', 1e-4);                                        %INPUT

  [timeGapClosure,y] = ode15s(@(t,x) funRHS(t,x), g.timeVector, y(end,:)', options);

  if timeGapClosure < g.timeEnd2
   % this is the gap closure... stop
     fprintf('Fuel-clad gap closed...\n');
     gap.open = 0;
     gap.clsd = 1;
     gap.depsz = clad.eps{3}(:,1) - fuel.epsT;
     gap.depsh = clad.eps{2}(:,1) - fuel.epsT;

   % Time vector for transient before gap closure
     g.timeStep = 0.1; %s
     g.timeVector = (timeGapClosure(end,1) :g.timeStep: g.timeEnd2);

     options = odeset('OutputFcn',@(t,y,flag) writeResults(t,y,flag), ...
                      'Mass', massMatrix, ...
                      'MaxStep', 1e-3, ...
                      'RelTol', 1e-6, ...                                  %INPUT
                      'AbsTol', 1e-4);                                     %INPUT

     ode15s(@(t,x) funRHS(t,x), g.timeVector, y(end,:)', options);
  end
  fprintf('Finished. Now run the script plotResults.m\n');
end