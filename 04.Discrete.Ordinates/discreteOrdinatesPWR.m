% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates the neutron transport in a 2D (x,y) unit cell 
% using the discrete ordinates method.

function discreteOrdinatesPWR

% Start stopwatch
  tic;

% Global structure with geometry parameters, total cross section and number
% of groups
  global g SigT ng
  
%--------------------------------------------------------------------------
% Path to macroscopic cross section data:
  path(path,['..' filesep '02.Macro.XS.421g']);
% Fill the structures fuel, clad and cool with the cross sections data
  fuel = macro421_UO2_03__900K;                                            % INPUT
  clad = macro421_Zry__600K;                                               % INPUT
  cool = macro421_H2OB__600K;                                              % INPUT

% Number of energy groups
  ng = fuel.ng;

%--------------------------------------------------------------------------
% Number of nodes
  g.nNodesX = 10;                                                          % INPUT
  g.nNodesY = 2;                                                           % INPUT

% Define the mesh step, nodes coordinates and node volumes
  g.delta = 0.2; %cm
  for iy = 1:g.nNodesY
      for ix = 1:g.nNodesX
          volume(ix,iy) = g.delta^2;
          if ix == 1 || ix == g.nNodesX
             volume(ix,iy) = volume(ix,iy) / 2;
          end
          if iy == 1 || iy == g.nNodesY
             volume(ix,iy) = volume(ix,iy) / 2;
          end
      end
  end

% define the material for each node (0 is coolant, 1 is cladding, 2 is fuel)
  mat = [2 2 2 2 2 1 0 0 0 0; ...
		 2 2 2 2 2 1 0 0 0 0 ]';

%--------------------------------------------------------------------------
% Path to Lebedev quadrature function:
  path(path,['..' filesep '00.Lebedev']);

% Number of discrete ordinates, an even integer (possible values are
% determined by the Lebedev quadratures: 6, 14, 26, 38, 50, 74, 86, 110,
% 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730,
% 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810):
  g.N = 110;                                                               % INPUT

% Get leb.x, leb.y and leb.z values for the g.N base points on a unit
% sphere as well as the associated weights leb.w (the unit sphere area
% corresponding to the base points, sum up to 4*pi) using the Lebedev
% quadrature rules.
  leb = getLebedevSphere(g.N);
  g.muX = leb.x;
  g.muY = leb.y;
  g.muZ = leb.z;
  g.W = leb.w;

% Find the reflective directions for X, Y and Z directions  
  for n = 1:g.N
      for nn = 1:g.N % find a symmetric direction
          if g.muX(nn) == -g.muX(n) && g.muY(nn) == g.muY(n) && g.muZ(nn) == g.muZ(n)
             g.nRefX(n) = nn;
          end
          if g.muY(nn) == -g.muY(n) && g.muX(nn) == g.muX(n) && g.muZ(nn) == g.muZ(n)
             g.nRefY(n) = nn;
          end
          if g.muZ(nn) == -g.muZ(n) && g.muX(nn) == g.muX(n) && g.muY(nn) == g.muY(n)
             g.nRefZ(n) = nn;
          end
      end
  end

%--------------------------------------------------------------------------
% Scattering source anisotropy: 0 -- P0 (isotropic), 1 -- P1
  g.L = 0;                                                                 % INPUT
% Calculate spherical harmonics for every ordinate
  for n = 1:g.N
      for jLgn = 0:g.L
          for m = -jLgn:+jLgn
              if jLgn == 0 && m == 0
                    g.R{n}(1+jLgn,(1+jLgn)+m) = 1;
              elseif jLgn == 1 && m == -1
                    g.R{n}(1+jLgn,(1+jLgn)+m) = g.muZ(n);
              elseif jLgn == 1 && m == 0
                    g.R{n}(1+jLgn,(1+jLgn)+m) = g.muX(n);
              elseif jLgn == 1 && m == 1
                    g.R{n}(1+jLgn,(1+jLgn)+m) = g.muY(n);
              end % of if
          end % of m
      end % of jLgn
  end % of cycle over ordinates
  
%--------------------------------------------------------------------------
% Construct the cross sections
  for iy = 1:g.nNodesY
      for ix = 1:g.nNodesX
          switch mat(ix,iy)
              case 2 % fuel
                   SigA{ix,iy} = fuel.SigF + fuel.SigC + fuel.SigL + sum(fuel.Sig2,2)';
                   for jLgn=0:g.L
                       SigS{1+jLgn}{ix,iy} = fuel.SigS{1+jLgn};
                   end
                   Sig2{ix,iy} = fuel.Sig2;
                   SigP{ix,iy} = fuel.SigP;
                   chi{ix,iy} = fuel.chi;
              case 1 % cladding
                   SigA{ix,iy} = clad.SigF + clad.SigC + clad.SigL + sum(clad.Sig2,2)';
                   for jLgn=0:g.L
                       SigS{1+jLgn}{ix,iy} = clad.SigS{1+jLgn};
                   end
                   Sig2{ix,iy} = clad.Sig2;
                   SigP{ix,iy} = clad.SigP;
                   chi{ix,iy} = clad.chi;
              case 0 % coolant
                   SigA{ix,iy} = cool.SigF + cool.SigC + cool.SigL + sum(cool.Sig2,2)';
                   for jLgn=0:g.L
                       SigS{1+jLgn}{ix,iy} = cool.SigS{1+jLgn};
                   end
                   Sig2{ix,iy} = cool.Sig2;
                   SigP{ix,iy} = cool.SigP;
                   chi{ix,iy} = cool.chi;
          end
          
        % Total cross section
          SigT{ix,iy} = SigA{ix,iy} + sum(SigS{1+0}{ix,iy},2)';
      end
  end

%--------------------------------------------------------------------------
% Count the number of equations.
  nEq = 0;
  for iy = 1:g.nNodesY
      for ix = 1:g.nNodesX
          for n = 1:g.N
              if g.muZ(n) >= 0  && ~(ix==1 && g.muX(n)>0) && ~(ix==g.nNodesX && g.muX(n)<0) ...
			                    && ~(iy==1 && g.muY(n)>0) && ~(iy==g.nNodesY && g.muY(n)<0)
                 nEq = nEq + ng;
              end
          end
      end
  end  
%--------------------------------------------------------------------------
  keff = [ ];
  residual = [ ];
% Number of outer iterations
  numIter = 200;                                                           % INPUT
% Set the initial flux equal 1.
  solution = ones(nEq,1);
% Main iteration loop
  for nIter = 1:numIter
    %-----------------------------------------------------------------------
    % Make a guess for the solution
    % Just take the solution from previous iteration as a guess
      guess = solution;
    %-----------------------------------------------------------------------
    % Convert 1D guess vector to the array of angular flux fi
      fi = convert(guess);
    %-----------------------------------------------------------------------
      for iy = 1:g.nNodesY
          for ix = 1:g.nNodesX
            % convert angular flux fi to the Legendre moments fiL
              for jLgn = 0:g.L
                  for m = -jLgn:+jLgn
                      SUM = zeros(ng,1);
                      for n = 1:g.N
                          SUM = SUM + fi(:,n,ix,iy) * g.R{n}(1+jLgn,(1+jLgn)+m) * g.W(n);
                      end
                      fiL{ix,iy}(:,1+jLgn,(1+jLgn)+m) = SUM;
                  end
              end
            % Scalar flux
              FI{ix,iy} = fiL{ix,iy}(:,1+0,1+0);
          end
      end
    %-----------------------------------------------------------------------
    % pRate is total neutron production rate
      pRate = 0;
    % aRate is total neutron absorption rate
      aRate = 0;
      for iy = 1:g.nNodesY
          for ix = 1:g.nNodesX
              pRate = pRate + ( SigP{ix,iy} + 2*sum(Sig2{ix,iy},2)' ) * FI{ix,iy} * volume(ix,iy);
              aRate = aRate + SigA{ix,iy} * FI{ix,iy} * volume(ix,iy);
          end
      end
    % We evaluate the multiplication factor as a ratio of neutron
    % production rate and neutron absorption rate (there is no neutron
    % leakage in the infinite lattice):
      keff = [keff pRate/aRate];
      fprintf('keff = %9.5f #OuterIter =%4i ',keff(end),nIter);

    %-----------------------------------------------------------------------
    % Calculate fission, (n,2n) and scattering neutron sources
      nEq = 0;
      for iy = 1:g.nNodesY
          for ix = 1:g.nNodesX
            % Fission source (1/s-cm3-steradian)
              qF = chi{ix,iy}' * SigP{ix,iy} * FI{ix,iy} / keff(end) / (4*pi);
            % Isotropic source from (n,2n) (1/s-cm3-steradian)
              q2 = 2*Sig2{ix,iy}' * FI{ix,iy} / (4*pi);
              for n = 1:g.N
                  if g.muZ(n) >= 0 && ~(ix==1 && g.muX(n)>0) && ~(ix==g.nNodesX && g.muX(n)<0) ...
				                   && ~(iy==1 && g.muY(n)>0) && ~(iy==g.nNodesY && g.muY(n)<0)
                   % Scattering source (1/s-cm3-steradian), isotropic (g.L =
                   % 0) or anisotropic (g.L > 0)
				     qS = zeros(ng,1);
                     for jLgn = 0:g.L
                         SUM = zeros(ng,1);
                         for m = -jLgn:+jLgn
                             SUM = SUM + fiL{ix,iy}(:,1+jLgn,(1+jLgn)+m) * g.R{n}(1+jLgn,(1+jLgn)+m);
                         end
                         qS = qS + (2*jLgn+1)*SigS{1+jLgn}{ix,iy}' * SUM / (4*pi);
                     end

                     nEq = nEq + 1;
                   % Right-hand side is a total neutron source:
                     qT(:,nEq) = qF + q2 + qS;
                  end
              end
          end
      end
      RHS = reshape(qT,[],1);

    %-----------------------------------------------------------------------
    % Relative residual reduction factor
      errtol = 1.e-4;                                                      % INPUT
    % maximum number of iterations
      maxit = 2000;                                                        % INPUT
    % Solver of a system of linear algebraic equations:
    % http://www4.ncsu.edu/~ctk/matlab_roots.html
      [solution, r, nInner] = bicgstab_(guess, RHS, @funDO, [errtol maxit]);

    % Save relative residual
      residual = [residual; r(end)/norm(RHS)];
      fprintf('numInner =%4i residual =%12.5e target =%12.5e\n',nInner, residual(end), errtol);
      
      if nInner == 0, break, end;

  end % of the main iteration loop
  
%--------------------------------------------------------------------------
% Find integral scalar flux in fuel, cladding and coolant (average spectra)
  FIFuel = zeros(ng,1);
  FIClad = zeros(ng,1);
  FICool = zeros(ng,1);
  for iy=1:g.nNodesY
      for ix=1:g.nNodesX
          switch mat(ix,iy)
              case 2
                    FIFuel = FIFuel + FI{ix,iy}*volume(ix,iy);
              case 1
                    FIClad = FIClad + FI{ix,iy}*volume(ix,iy);
              case 0
                    FICool = FICool + FI{ix,iy}*volume(ix,iy);
          end
      end
  end
%--------------------------------------------------------------------------
% Write down the results  
  fileID = fopen('resultsPWR.m','w');
     fprintf(fileID,'function s = resultsPWR\n\n');
     fprintf(fileID,'%% Results for 2D neutron transport calculation in the PWR-like unit cell using method of discrete ordinates\n');
     fprintf(fileID,'%% Number of ordinates used is %3i\n', g.N);
     fprintf(fileID,'%% Scattering source anisotropy approximation is P%i\n\n', g.L);
   % Stop stopwatch
     elapsedTime = toc;
     fprintf(fileID,'%% elapsedTime = %12.5e s / %12.5e min / %12.5e hrs;\n\n',elapsedTime, elapsedTime/60, elapsedTime/3600);
     fprintf(fileID,'s.keff = %12.5e;\n\n',keff(end));
     fprintf(fileID,'s.keffHistory = ['); fprintf(fileID,'%12.5e ',keff); fprintf(fileID,'];\n\n');
     fprintf(fileID,'s.residualHistory = ['); fprintf(fileID,'%12.5e ',residual); fprintf(fileID,'];\n\n');

     fprintf(fileID,'s.x = ['); fprintf(fileID,'%e ',(0:g.delta:g.delta*(g.nNodesX-1))); fprintf(fileID,'];\n\n');% mesh (cm)
     fprintf(fileID,'s.eg = ['); fprintf(fileID,'%e ',(fuel.eg(1:ng)+fuel.eg(2:ng+1))/2.0); fprintf(fileID,'];\n\n');% group energies
     du = log( fuel.eg(2:ng+1)./fuel.eg(1:ng) )';
     fprintf(fileID,'%% Neutron spectrum in fuel, cladding and coolant:\n');
     fprintf(fileID,'s.FIFuel_du = [');  fprintf(fileID,'%e ',FIFuel(1:ng) ./ du);  fprintf(fileID,'];\n');
     fprintf(fileID,'s.FIClad_du = [');  fprintf(fileID,'%e ',FIClad(1:ng) ./ du);  fprintf(fileID,'];\n');
     fprintf(fileID,'s.FICool_du = [');  fprintf(fileID,'%e ',FICool(1:ng) ./ du);  fprintf(fileID,'];\n\n');
     
     fprintf(fileID,'%% Thermal-, resonance- and fast-group fluxes along the cell centerline:\n');
     fprintf(fileID,'s.FI_T = [');
     for ix = 1:g.nNodesX
         fprintf(fileID,'%e ',sum(FI{ix,1}(1:50)));  % < 1 eV
     end
     fprintf(fileID,'];\n');
     fprintf(fileID,'s.FI_R = [');
     for ix = 1:g.nNodesX
         fprintf(fileID,'%e ',sum(FI{ix,1}(51:287)));  % < 0.1 Mev
     end
     fprintf(fileID,'];\n');
     fprintf(fileID,'s.FI_F = [');
     for ix = 1:g.nNodesX
         fprintf(fileID,'%e ',sum(FI{ix,1}(288:421)));  % > 0.1 MeV
     end
     fprintf(fileID,'];\n\n');
     fprintf(fileID,'end');
  fclose(fileID);
  
%--------------------------------------------------------------------------
% Plot the mesh
  plot2D(g.nNodesX, g.nNodesY, g.delta, mat, 'Unit cell: materials', 'Fig_01_discreteOrdinates.pdf');

%--------------------------------------------------------------------------
% Plot the results
  s = resultsPWR;

  f = figure('visible','off');
  plot(keff,'-or');
  ylim(ylim);
  grid on;
  xlabel('Iteration number');
  ylabel('k-effective');
  saveas(f, 'Fig_02_discreteOrdinates.pdf');

  f = figure('visible','off');
  semilogy(residual,'-or');
  ylim(ylim);
  grid on;
  xlabel('Iteration number');
  ylabel('Relative residual error');
  saveas(f, 'Fig_03_discreteOrdinates.pdf');
  
  f = figure('visible','off');
  semilogx(s.eg,s.FIFuel_du,'-r',...
           s.eg,s.FIClad_du,'-g',...
           s.eg,s.FICool_du,'-b');
  grid on;
  xlabel('Energy (eV)');
  ylabel('Neutron flux per unit lethargy (a.u.)');
  legend('Fuel','Cladding','Coolant','Location','northwest');
  saveas(f, 'Fig_04_discreteOrdinates.pdf');
  
  f = figure('visible','off');
  plot(s.x,s.FI_F,'-or',...
       s.x,s.FI_R,'-og', ...
       s.x,s.FI_T,'-ob');
  ylim(ylim);
  grid on;
  xlabel('Distance from the cell centre (cm)');
  ylabel('Neutron flux (a.u.)');
  legend('Fast','Resonance','Thermal','Location','northwest');
  saveas(f, 'Fig_05_discreteOrdinates.pdf');
  
  for iy=1:g.nNodesY
      for ix=1:g.nNodesX
          funT(ix,iy) = sum(FI{ix,iy}(1:50),1);
          funR(ix,iy) = sum(FI{ix,iy}(51:355),1);
          funF(ix,iy) = sum(FI{ix,iy}(356:421),1);
      end
  end
  plot2D(g.nNodesX, g.nNodesY, g.delta, funT, 'Thermal flux distribution', 'Fig_06_discreteOrdinates.pdf');
  plot2D(g.nNodesX, g.nNodesY, g.delta, funR, 'Resonance flux distribution', 'Fig_07_discreteOrdinates.pdf');
  plot2D(g.nNodesX, g.nNodesY, g.delta, funF, 'Fast flux distribution', 'Fig_08_discreteOrdinates.pdf');
    
end