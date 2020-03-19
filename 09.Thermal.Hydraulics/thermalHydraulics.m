% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates coupled thermal hydraulics and fuel behaviour for
% one PWR-like subassembly under 1D approximation under LOCA conditions

function thermalHydraulics

  clear all;

% Global structures
  global g fuel clad cool

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
  initializeCoolant

%--------------------------------------------------------------------------
% Time vector
  g.timeStep = 1; %s
  g.timeEnd = 600;  %s

  g.timeVector = (0 :g.timeStep: g.timeEnd);
  
%--------------------------------------------------------------------------
% Diagonal "mass matrix" to multiply the time derivatives of the ODE
% system. Needed by ode15s to solve together ODEs (1) and AEs (0) 
  massMatrix = diag([... 
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
                   'Events',@cladFailureEvent, ...
                   'RelTol', 1e-6, ...                                     %INPUT
                   'AbsTol', 1e-4);                                        %INPUT
  
%--------------------------------------------------------------------------
% Initial guess for unknowns
  y = [...
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
  [g.timeEnd1,y] = ode15s(@(t,x) funRHS(t,x), g.timeVector, y, options);

  if g.timeEnd1(end) < g.timeEnd
   % this is the clad failure... go on with the failed clad
     fprintf('clad failure at time %9.2f\n',g.timeEnd1(end));
   % clad failure time
     clad.failTime = g.timeEnd1(end);
   % clad failure index
     clad.fail = true;
     options = odeset('OutputFcn',@(t,y,flag) writeResults(t,y,flag), ...
                      'Mass', massMatrix, ...
                      'MaxStep', 1e-3, ...
                      'RelTol', 1e-6, ...                                  %INPUT
                      'AbsTol', 1e-4);                                     %INPUT

     ode15s(@(t,x) funRHS(t,x), (clad.failTime :g.timeStep: g.timeEnd), y(end,:)', options);
  end
  fprintf('Finished. Now run the script plotResults.m\n');
end