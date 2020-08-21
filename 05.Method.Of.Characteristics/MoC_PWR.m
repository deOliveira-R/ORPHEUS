% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates the neutron transport in a 2D (x,y) unit cell 
% using the method of characteristics.

function MoC_PWR

% Start stopwatch
  tic;

% Global structure with geometry parameters g, flux fi, total cross section
% SigT, neutron source and direction IDs
  global g fi SigT q d

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
% Number of rays (must be 8 according to the selected mesh)
  nRays = 8;
% direction IDs
  d.EAST=1; d.NORTH_EAST=2; d.NORTH=3; d.NORTH_WEST=4; d.WEST=5; d.SOUTH_WEST=6; d.SOUTH=7; d.SOUTH_EAST=8;

%--------------------------------------------------------------------------  
% Number of nodes
  g.nNodes = 10;                                                           % INPUT

% Define the mesh step, nodes coordinates and node volumes
  g.delta = 0.2; %cm                                                       % INPUT
  for iy = 1:g.nNodes
      for ix = 1:g.nNodes
          volume(ix,iy) = g.delta^2;
          if ix == 1 || ix == g.nNodes
             volume(ix,iy) = volume(ix,iy) / 2;
          end
          if iy == 1 || iy == g.nNodes
             volume(ix,iy) = volume(ix,iy) / 2;
          end
      end
  end

% define the material for each node (0 is coolant, 1 is cladding, 2 is fuel)
  mat = [2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0; ...
         2 2 2 2 2 1 0 0 0 0 ]';

%--------------------------------------------------------------------------
% Construct the cross sections and set the initial flux to 1
  for iy = 1:g.nNodes
      for ix = 1:g.nNodes
          switch mat(ix,iy)
              case 2 % fuel
                   SigA{ix,iy} = fuel.SigF + fuel.SigC + fuel.SigL + sum(fuel.Sig2,2)';
                   SigS0{ix,iy} = fuel.SigS{1+0};
                   Sig2{ix,iy} = fuel.Sig2;
                   SigP{ix,iy} = fuel.SigP;
                   chi{ix,iy} = fuel.chi;
              case 1 % cladding
                   SigA{ix,iy} = clad.SigF + clad.SigC + clad.SigL + sum(clad.Sig2,2)';
                   SigS0{ix,iy} = clad.SigS{1+0};
                   Sig2{ix,iy} = clad.Sig2;
                   SigP{ix,iy} = clad.SigP;
                   chi{ix,iy} = clad.chi;
              case 0 % coolant
                   SigA{ix,iy} = cool.SigF + cool.SigC + cool.SigL + sum(cool.Sig2,2)';
                   SigS0{ix,iy} = cool.SigS{1+0};
                   Sig2{ix,iy} = cool.Sig2;
                   SigP{ix,iy} = cool.SigP;
                   chi{ix,iy} = cool.chi;
          end

        % Total cross section
          SigT{ix,iy} = SigA{ix,iy} + sum(SigS0{ix,iy},2)';

        % Set the initial values for nRays components of the angular flux:
          fi{ix,iy} = ones(ng,nRays);
      end
  end
  

%--------------------------------------------------------------------------
  keff = [ ];
% Number of outer iterations
  numIter = 200;                                                           % INPUT
% Main iteration loop
  for nIter = 1:numIter
    %-----------------------------------------------------------------------
    % pRate is total neutron production rate
      pRate = 0;
    % aRate is total neutron absorption rate
      aRate = 0;
      for iy = 1: g.nNodes
          for ix = 1:g.nNodes
            % Scalar flux
              FI{ix,iy} = sum(fi{ix,iy},2);
              pRate = pRate + ( SigP{ix,iy} + 2*sum(Sig2{ix,iy},2)' ) * FI{ix,iy} * volume(ix,iy);
              aRate = aRate + SigA{ix,iy} * FI{ix,iy} * volume(ix,iy);
          end
      end
    % We evaluate the multiplication factor as a ratio of neutron
    % production rate and neutron absorption rate (there is no neutron
    % leakage in the infinite lattice):
      keff = [keff pRate/aRate];
      fprintf('keff = %9.5f #nOuter = %3i\n',keff(end),nIter);
      
    %----------------------------------------------------------------------
    % Calculate the distribution of the isotropic neutron source per ray
      for iy = 1:g.nNodes
          for ix = 1:g.nNodes
              q{ix,iy} = (chi{ix,iy}'*SigP{ix,iy}/keff(end) + 2*Sig2{ix,iy}' + SigS0{ix,iy}') * FI{ix,iy};
              q{ix,iy} = q{ix,iy} / nRays;
          end
      end
    %----------------------------------------------------------------------
    % Go along the horizontal rays east, reflect, return west, reflect
      for iy = 1:g.nNodes
          for ix = 1:g.nNodes
              flyFrom(ix, iy, d.EAST);
          end
          fi{g.nNodes,iy}(:,d.WEST) = fi{g.nNodes,iy}(:,d.EAST); % reflect on the right boundary
          for ix = g.nNodes:-1:1
              flyFrom(ix, iy, d.WEST);
          end
          fi{1,iy}(:,d.EAST) = fi{1,iy}(:,d.WEST); % reflect on the left boundary
      end
    %----------------------------------------------------------------------
    % Go along the vertical rays south, reflect, return north, reflect
      for ix = 1:g.nNodes
          for iy = 1:g.nNodes
              flyFrom(ix, iy, d.SOUTH);
          end
          fi{ix,g.nNodes}(:,d.NORTH) = fi{ix,g.nNodes}(:,d.SOUTH); % reflect on the bottom boundary
          for iy = g.nNodes:-1:1
              flyFrom(ix, iy, d.NORTH);
          end
          fi{ix,1}(:,d.SOUTH) = fi{ix,1}(:,d.NORTH); % reflect on the bottom boundary
      end
    %----------------------------------------------------------------------
    % Go along the diagonal ray south-east, reflect, return north-west, reflect
      ixx = 1;
      iyy = 1;
      while ixx < g.nNodes
          flyFrom(ixx, iyy, d.SOUTH_EAST);
          ixx = ixx + 1;
          iyy = iyy + 1;
      end
      fi{g.nNodes,g.nNodes}(:,d.NORTH_WEST) = fi{g.nNodes,g.nNodes}(:,d.SOUTH_EAST); % reflect on the bottom right corner
      fi{g.nNodes,g.nNodes}(:,d.NORTH_EAST) = fi{g.nNodes,g.nNodes}(:,d.SOUTH_EAST); % reflect on the bottom right corner
      fi{g.nNodes,g.nNodes}(:,d.SOUTH_WEST) = fi{g.nNodes,g.nNodes}(:,d.SOUTH_EAST); % reflect on the bottom right corner
      ixx = g.nNodes;
      iyy = g.nNodes;
      while ixx > 1
          flyFrom(ixx, iyy, d.NORTH_WEST);
          ixx = ixx - 1;
          iyy = iyy - 1;
      end
      fi{1,1}(:,d.SOUTH_EAST) = fi{1,1}(:,d.NORTH_WEST); % reflect on the top left corner
      fi{1,1}(:,d.NORTH_EAST) = fi{1,1}(:,d.NORTH_WEST); % reflect on the top left corner
      fi{1,1}(:,d.SOUTH_WEST) = fi{1,1}(:,d.NORTH_WEST); % reflect on the top left corner
    %----------------------------------------------------------------------
    % Go along the diagonal ray south-west, reflect and return north-east
      ixx = g.nNodes;
      iyy = 1;
      while ixx > 1 
          flyFrom(ixx, iyy, d.SOUTH_WEST);
          ixx = ixx - 1;
          iyy = iyy + 1;
      end
      fi{1,g.nNodes}(:,d.NORTH_EAST) = fi{1,g.nNodes}(:,d.SOUTH_WEST); % reflect on the bottom left corner
      fi{1,g.nNodes}(:,d.NORTH_WEST) = fi{1,g.nNodes}(:,d.SOUTH_WEST); % reflect on the bottom left corner
      fi{1,g.nNodes}(:,d.SOUTH_EAST) = fi{1,g.nNodes}(:,d.SOUTH_WEST); % reflect on the bottom left corner
      ixx = 1;
      iyy = g.nNodes;
      while ixx < g.nNodes 
          flyFrom(ixx, iyy, d.NORTH_EAST);
          ixx = ixx + 1;
          iyy = iyy - 1;
      end
      fi{g.nNodes,1}(:,d.SOUTH_WEST) = fi{g.nNodes,1}(:,d.NORTH_EAST); % reflect on the top right corner
      fi{g.nNodes,1}(:,d.SOUTH_EAST) = fi{g.nNodes,1}(:,d.NORTH_EAST); % reflect on the top right corner
      fi{g.nNodes,1}(:,d.NORTH_WEST) = fi{g.nNodes,1}(:,d.NORTH_EAST); % reflect on the top right corner
    %----------------------------------------------------------------------
    % Go along non-diagonal rays south-east rays and reflect on the right boundary
      for ix=g.nNodes-1:-1:2
          ixx = ix;
          iyy = 1;
          while ixx < g.nNodes
                flyFrom(ixx, iyy, d.SOUTH_EAST);
                ixx = ixx + 1;
                iyy = iyy + 1;
          end
          fi{ixx,iyy}(:,d.SOUTH_WEST) = fi{ixx,iyy}(:,d.SOUTH_EAST); % reflect on the right boundary
      end
      
    % Go along non-diagonal rays south-east rays and reflect on the bottom boundary
      for iy=2:g.nNodes-1
          ixx = 1;
          iyy = iy;
          while iyy < g.nNodes
              flyFrom(ixx, iyy, d.SOUTH_EAST);
              ixx = ixx + 1;
              iyy = iyy + 1;
          end
          fi{ixx,iyy}(:,d.NORTH_EAST) = fi{ixx,iyy}(:,d.SOUTH_EAST); % reflect on the bottom boundary
      end
    %----------------------------------------------------------------------
    % Go along non-diagonal rays south-west rays and reflect on the left boundary
      for ix = 2:g.nNodes-1
          ixx = ix;
          iyy = 1;
          while ixx > 1
                flyFrom(ixx, iyy, d.SOUTH_WEST);
                ixx = ixx - 1;
                iyy = iyy + 1;
          end
          fi{ixx,iyy}(:,d.SOUTH_EAST) = fi{ixx,iyy}(:,d.SOUTH_WEST); % reflect on the left boundary
      end
      
    % Go along non-diagonal rays south-west rays and reflect on the bottom boundary
      for iy = 2: g.nNodes-1
          ixx = g.nNodes;
          iyy = iy;
          while iyy < g.nNodes
                flyFrom(ixx, iyy, d.SOUTH_WEST);
                ixx = ixx - 1;
                iyy = iyy + 1;
          end
          fi{ixx,iyy}(:,d.NORTH_WEST) = fi{ixx,iyy}(:,d.SOUTH_WEST); % reflect on the bottom boundary
      end
    %----------------------------------------------------------------------
    % Go along non-diagonal rays north-west and reflect on the left boundary
      for ix = 2:g.nNodes-1
          ixx = ix;
          iyy = g.nNodes;
          while ixx > 1
                flyFrom(ixx, iyy, d.NORTH_WEST);
                ixx = ixx - 1;
                iyy = iyy - 1;
          end
          fi{ixx,iyy}(:,d.NORTH_EAST) = fi{ixx,iyy}(:,d.NORTH_WEST); % reflect on the left boundary
      end
      
    % Go along non-diagonal rays north-west and reflect on the top boundary
      for iy = g.nNodes-1:-1:2
          ixx = g.nNodes;
          iyy = iy;
          while iyy > 1
                flyFrom(ixx, iyy, d.NORTH_WEST);
                ixx = ixx - 1;
                iyy = iyy - 1;
          end
          fi{ixx,iyy}(:,d.SOUTH_WEST) = fi{ixx,iyy}(:,d.NORTH_WEST); % reflect on the top boundary
      end
    %----------------------------------------------------------------------
    % Go along non-diagonal rays north-east and reflect on the top boundary
      for iy = 2:g.nNodes-1
          ixx = 1;
          iyy = iy;
          while iyy > 1
                flyFrom(ixx, iyy, d.NORTH_EAST);
                ixx = ixx + 1;
                iyy = iyy - 1;
          end
          fi{ixx,iyy}(:,d.SOUTH_EAST) = fi{ixx,iyy}(:,d.NORTH_EAST); % reflect on the top boundary
      end
      
    % Go along non-diagonal rays north-east and reflect on the right boundary
      for ix = 2:g.nNodes-1
          ixx = ix;
          iyy = g.nNodes;
          while ixx < g.nNodes
                flyFrom(ixx, iyy, d.NORTH_EAST);
                ixx = ixx + 1;
                iyy = iyy - 1;
          end
          fi{ixx,iyy}(:,d.NORTH_WEST) = fi{ixx,iyy}(:,d.NORTH_EAST); % reflect on the right boundary
      end

    %----------------------------------------------------------------------
    % Normalize the flux so that the total production rate is 100 n/s
      for iy=1:g.nNodes
          for ix=1:g.nNodes
              fi{ix,iy} = fi{ix,iy} / pRate * 100;
          end
      end

  end % of the main iteration loop
%--------------------------------------------------------------------------
% Find integral scalar flux in fuel, cladding and coolant (average spectra)
  vol_fuel = sum(volume(1:5,:),'all');
  vol_clad = sum(volume(6,:),'all');
  vol_cool = sum(volume(7:end,:),'all');

  FIFuel = zeros(ng,1);
  FIClad = zeros(ng,1);
  FICool = zeros(ng,1);
  for iy=1:g.nNodes
      for ix=1:g.nNodes
          switch mat(ix,iy)
              case 2
                    FIFuel = FIFuel + FI{ix,iy}*volume(ix,iy)/vol_fuel;
              case 1
                    FIClad = FIClad + FI{ix,iy}*volume(ix,iy)/vol_clad;
              case 0
                    FICool = FICool + FI{ix,iy}*volume(ix,iy)/vol_cool;
          end
      end
  end
%--------------------------------------------------------------------------
% Write down the results  
  fileID = fopen('resultsPWR.m','w');
     fprintf(fileID,'function s = resultsPWR\n\n');
     fprintf(fileID,'%% Results for 2D neutron transport calculation in PWR-like unit cell using method of characteristics\n\n');     
   % Stop stopwatch
     elapsedTime = toc;
     fprintf(fileID,'%% elapsedTime = %12.5e s / %12.5e min / %12.5e hrs;\n\n',elapsedTime, elapsedTime/60, elapsedTime/3600);
     fprintf(fileID,'s.keff = %12.5e;\n\n',keff(end));
     fprintf(fileID,'s.keffHistory = ['); fprintf(fileID,'%12.5e ',keff); fprintf(fileID,'];\n\n');     
     
     fprintf(fileID,'s.x = ['); fprintf(fileID,'%e ',(0:g.delta:g.delta*(g.nNodes-1))); fprintf(fileID,'];\n\n');% mesh (cm)
     fprintf(fileID,'s.eg = ['); fprintf(fileID,'%e ',(fuel.eg(1:ng)+fuel.eg(2:ng+1))/2.0); fprintf(fileID,'];\n\n');% group energies
     du = log( fuel.eg(2:ng+1)./fuel.eg(1:ng) )';
     fprintf(fileID,'%% Neutron spectrum in fuel, cladding and coolant:\n');     
     fprintf(fileID,'s.FIFuel_du = [');  fprintf(fileID,'%e ',FIFuel(1:ng) ./ du);  fprintf(fileID,'];\n');
     fprintf(fileID,'s.FIClad_du = [');  fprintf(fileID,'%e ',FIClad(1:ng) ./ du);  fprintf(fileID,'];\n');
     fprintf(fileID,'s.FICool_du = [');  fprintf(fileID,'%e ',FICool(1:ng) ./ du);  fprintf(fileID,'];\n\n');
     
     fprintf(fileID,'%% Thermal-, resonance- and fast-group fluxes along the cell centerline:\n');
     fprintf(fileID,'s.FI_T = [');
     for ix = 1:g.nNodes
         fprintf(fileID,'%e ',sum(FI{ix,1}(1:50)));  % < 1 eV
     end
     fprintf(fileID,'];\n');
     fprintf(fileID,'s.FI_R = [');
     for ix = 1:g.nNodes
         fprintf(fileID,'%e ',sum(FI{ix,1}(51:287)));  % < 0.1 Mev
     end
     fprintf(fileID,'];\n');
     fprintf(fileID,'s.FI_F = [');
     for ix = 1:g.nNodes
         fprintf(fileID,'%e ',sum(FI{ix,1}(288:421))); % > 0.1 MeV
     end
     fprintf(fileID,'];\n\n');
     fprintf(fileID,'end');
  fclose(fileID);
  
%--------------------------------------------------------------------------
% Plot the rays
  plotRays(g.nNodes, g.delta, 'Rays used to track neutrons', 'MOC_01_tracks.pdf');
% Plot the materials
  plot2D(g.nNodes, g.delta, mat, 'Materials', 'MOC_02_mesh.pdf');

%--------------------------------------------------------------------------
% Plot the results
  s = resultsPWR;
  
  f = figure('visible','off');
  plot(keff,'-or');
  ylim(ylim);
  grid on;
  xlabel('Iteration number');
  ylabel('k-effective');
  saveas(f, 'MOC_03_keff.pdf');
  
  f = figure('visible','off');
  semilogx(s.eg,s.FIFuel_du,'-r',...
           s.eg,s.FIClad_du,'-g',...
           s.eg,s.FICool_du,'-b')
  grid on;
  xlabel('Energy (eV)');
  ylabel('Neutron flux per unit lethargy (a.u.)');
  legend('Fuel','Cladding','Coolant','Location','best')
  saveas(f, 'MOC_04_flux_lethargy.pdf');
  
  f = figure('visible','off');
  plot(s.x,s.FI_F,'-or',...
       s.x,s.FI_R,'-og',...
       s.x,s.FI_T,'-ob')
  ylim(ylim);
  grid on;
  xlabel('Distance from the cell centre (cm)');
  ylabel('Neutron flux (a.u.)');
  legend('Fast','Resonance','Thermal','Location','northwest')
  saveas(f, 'MOC_05_flux_cell.pdf');
  
  for ix=1:g.nNodes
      for iy=1:g.nNodes
          funT(ix,iy) = sum(FI{ix,iy}(1:50),1);
          funR(ix,iy) = sum(FI{ix,iy}(51:355),1);
          funF(ix,iy) = sum(FI{ix,iy}(356:421),1);
      end
  end
  plot2D(g.nNodes, g.delta, funT, 'Thermal flux distribution', 'MOC_06_flux_thermal.pdf');
  plot2D(g.nNodes, g.delta, funR, 'Resonance flux distribution', 'MOC_07_flux_resonance.pdf');
  plot2D(g.nNodes, g.delta, funF, 'Fast flux distribution', 'MOC_08_flux_fast.pdf');

end
