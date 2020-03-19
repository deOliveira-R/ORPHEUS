% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates the neutron diffusion in a 1D PWR-like
% subassembly

function CORE1D

  global nReflBot nFuel nReflTop nNodes nEdges dz Sig D_ nGroups;
  
%--------------------------------------------------------------------------
% Number of energy groups
%--------------------------------------------------------------------------
% first group is fast; second group is thermal
  nGroups = 2;

%--------------------------------------------------------------------------
% Geometry specification
%--------------------------------------------------------------------------
% Approximate PWR square subassembly dimensions:
% bottom axial reflector height
  botReflHeight = 50; %cm                                            %INPUT
% fuel height
  fuelHeight = 300; %cm                                              %INPUT
% top axial reflector height
  topReflHeight = 50; %cm                                            %INPUT
% subassembly flat-to-flat distance
  f2f = 21.61; %cm                                                   %INPUT

% axial node height
  dz = 5; %cm                                                        %INPUT

% number of nodes  
  nReflBot = botReflHeight / dz; 
  nFuel = fuelHeight / dz; 
  nReflTop = topReflHeight / dz;  
  nNodes = nReflBot + nFuel + nReflTop; 
  nEdges = nNodes + 1;
  
% cross-sectional area of one SA (cm2)
  az = f2f^2;

% volume of one node (cm3):
  dv = az * dz;
  
% Average power of one fuel SA (W) is total PWR power divided by the number
% of fuel SAs:
  Power = 1.76752E+07;                                               %INPUT
  
%--------------------------------------------------------------------------
% Cross sections specification
%--------------------------------------------------------------------------
% Material 1: axial reflector
  axialRefl.Transport  = [0.3416;0.9431];                            %INPUT
  axialRefl.Absorption = [0.0029;0.0933];                            %INPUT
  axialRefl.Fission    = [0.0;0.0];                                  %INPUT
  axialRefl.Production = [0.0;0.0];                                  %INPUT
  axialRefl.Chi        = [0.0;0.0];                                  %INPUT
  axialRefl.Scattering = [2.4673e-04;0.0];                           %INPUT

% Material 2: fuel (3.1 w/o)
  fuel.Transport  = [0.2181;0.7850];                                 %INPUT
  fuel.Absorption = [0.0096;0.0959];                                 %INPUT
  fuel.Fission    = [0.0024;0.0489];                                 %INPUT
  fuel.Production = [0.0061;0.1211];                                 %INPUT
  fuel.Chi        = [1.0;0.0];                                       %INPUT
  fuel.Scattering = [0.0160;0.0];                                    %INPUT

% transport XS  
  Sig.T = [repmat(axialRefl.Transport,1,nReflBot),...
           repmat(fuel.Transport,1,nFuel),...
           repmat(axialRefl.Transport,1,nReflTop)];
% absorption XS  
  Sig.A = [repmat(axialRefl.Absorption,1,nReflBot),...
           repmat(fuel.Absorption,1,nFuel),...
           repmat(axialRefl.Absorption,1,nReflTop)];
% fission XS
  Sig.F = [repmat(axialRefl.Fission,1,nReflBot),...
           repmat(fuel.Fission,1,nFuel),...
           repmat(axialRefl.Fission,1,nReflTop)];
% production XS
  Sig.P = [repmat(axialRefl.Production,1,nReflBot),...
           repmat(fuel.Production,1,nFuel),...
           repmat(axialRefl.Production,1,nReflTop)];
% fission spectrum
  chi = [repmat(axialRefl.Chi,1,nReflBot),...
         repmat(fuel.Chi,1,nFuel),...
         repmat(axialRefl.Chi,1,nReflTop)];
% scattering XS
  Sig.S = [repmat(axialRefl.Scattering,1,nReflBot),...
           repmat(fuel.Scattering,1,nFuel),...
           repmat(axialRefl.Scattering,1,nReflTop)];
   
  for ig = 1:nGroups
      Sig.T_(ig,:) = [Sig.T(ig,1); interp1(Sig.T(ig,:),1.5:nNodes-0.5)'; Sig.T(ig,nNodes)];
  end

  D_ = 1 ./ (3*Sig.T_);
  
% Energy (J) released per fission
  ePerFission = [0.3213e-10; 0.3206e-10];
  
% Initialize scalar flux as a 2D matrix
  fi=ones(nGroups,nNodes);

% Initialize net current as a 2D matrix
  J=zeros(nGroups,nEdges);

% Make 1D solution vector
  solution = reshape(fi,[],1);
  
% inner and outer iteration counters
  nInner = 1;
  nOuter = 1;
  
% Loop for outer (power) iterations over the eigenvalue (k-eff)
  while nInner > 0

    % Node-wise production rate (1D vectors of length nNodes):
      pRate = sum(Sig.P .* fi * dv, 1);
    % Node-wise absorption rate (1D vector of length nNodes):
      aRate = sum(Sig.A .* fi * dv, 1);

    % Evaluate the gradient by a simple finite difference scheme
      dfidz = [(fi(:,1)-0)/(0.5*dz), diff(fi,1,2)/dz, (0-fi(:,nNodes))/(0.5*dz)];
    % Fick's law 
      J = -D_ .* dfidz;
    % Leakage rate (scalar):
      lRate = sum(-J(:,1)*az + J(:,nEdges)*az, 1);
          
    % We evaluate the multiplication factor as a ratio of neutron
    % production rate and neutron absorption + leakage rates
      keff = sum(pRate,2) / (sum(aRate,2) + lRate);

    %-----------------------------------------------------------------------
    % Make a guess for the solution
    % Just take the solution from previous iteration as a guess
      guess = solution;

    % A fission neutron source is a 1D column-vector [nGroups*nNodes,1]
      RHS = reshape(chi .* repmat(pRate/keff,nGroups,1) / keff, [], 1);

    %-----------------------------------------------------------------------
    % Relative residual reduction factor
      errtol = 1.e-6;                                                %INPUT
    % maximum number of iterations
      maxit = 2000;                                                  %INPUT
    % Solver of a system of linear algebraic equations:
    % http://www4.ncsu.edu/~ctk/matlab_roots.html
      [solution, r, nInner] = bicgstab_(guess, RHS, @funCORE1D, [errtol maxit]);
      
      nOuter = nOuter + 1;

      fi = reshape(solution,nGroups,nNodes);
      
    % Node-wise fission rate, 2D vectors (nGroups, nNodes):
      fRate = Sig.F .* fi * dv;
      pow = sum(sum(fRate,2) .* ePerFission,1);
    % Normalize the flux to the given Power
      fi = fi * Power/pow;
    
    % Print on the screen multiplication factor and outer iteration number
      fprintf('keff = %9.5f nOuter=%4i nInner=%4i residual=%12.5e\n',keff,nOuter,nInner,r(end)/norm(RHS));
        
  end
  
% Write results (scalar flux) to file results.m. To make plots run function
% makePlots.
  fdres=fopen('results.m','w');
  fprintf(fdres,'function s=results\n\n');
  fprintf(fdres,'s.keff=%e;\n\n',keff);
  
  %coordinates for flux
  zNodes=(dz/2 : dz: dz/2+dz*(nNodes-1));
  fprintf(fdres,'s.zNodes=[');
      fprintf(fdres,'%e ',zNodes);
  fprintf(fdres,'];\n\n');

  %coordinates for current
  zEdges=(0 : dz : dz*(nEdges-1));
  fprintf(fdres,'s.zEdges=[');
  fprintf(fdres,'%e ',zEdges);
  fprintf(fdres,'];\n\n');
  
  % flux
  fprintf(fdres,'s.fi=[...\n');
  for i=1:nGroups
      fprintf(fdres,'%e ',fi(i,1:nNodes));
      fprintf(fdres,'; ...\n');
  end
  fprintf(fdres,'];\n');
  
% net current
  fprintf(fdres,'s.J=[...\n');
  for i=1:nGroups
      fprintf(fdres,'%e ',J(i,1:nEdges));
      fprintf(fdres,'; ...\n');
  end
  fprintf(fdres,'];\n\n');
  
  fprintf(fdres,'end');
  
  %close results.m
  fclose(fdres);

%--------------------------------------------------------------------------
% Plot the results
  
  f = figure('visible','off');
  plot(zNodes,fi(1,:), '-or', ...
       zNodes,fi(2,:),'-ob');
  grid on;
  xlabel('z (cm)')
  ylabel('Neutron flux (n/cm2s)')
  legend('Fast','Thermal','Location','best')
  saveas(f, 'Fig_01_Neutron_Flux.pdf');
  
  f = figure('visible','off');
  plot(zEdges,J(1,:), '-or', ...
       zEdges,J(2,:),'-ob');
  grid on;
  xlabel('x (cm)')
  ylabel('Neutron net current (n/cm2s)')
  legend('Fast','Thermal','Location','best')
  saveas(f, 'Fig_02_Neutron_Net_Current.pdf');
end