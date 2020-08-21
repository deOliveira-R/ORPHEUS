% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function calculates the neutron transport in a 2D (x,y) unit cell
% similar to the unit cell of the pressurized water reactor using the Monte
% Carlo method.
%
% Credits: Inspired by the matlab functions mctptr and mctrk published in
% Applied Reactor Physics, Presses Internationales Polytechnique by Alain
% Hebert, Ecole Polytechnique de Montreal, based on a Fortran program
% written by Benoit Arsenault.

function monteCarloPWR

  clear all;

% Start stopwatch
  tic;

%--------------------------------------------------------------------------
% Number of source neutrons
  numNeutrons_born = 100;                                                  % INPUT

% Number of inactive source cycles to skip before starting k-eff
% accumulation
  numCycles_inactive = 100;                                                % INPUT
  
% Number of active source cycles for k-eff accumulation
  numCycles_active = 2000;                                                 % INPUT

% Size of the square unit cell
  pitch = 3.6; %cm                                                         % INPUT

%--------------------------------------------------------------------------
% Path to macroscopic cross section data:
  path(path,['..' filesep '02.Macro.XS.421g']);
% Fill the structures fuel, clad and cool with the cross sections data
  fuel = macro421_UO2_03__900K;                                            % INPUT
  clad = macro421_Zry__600K;                                               % INPUT
  cool = macro421_H2OB__600K;                                              % INPUT

% Define the majorant: the maximum total cross section vector 
  SigTmax = max([fuel.SigT; clad.SigT; cool.SigT]);

% Number of energy groups
  ng = fuel.ng;

%--------------------------------------------------------------------------
% Detectors
  detectS = zeros(1,ng);

%--------------------------------------------------------------------------
% Four main vectors describing the neutrons in a batch 
  x = zeros(1,numNeutrons_born*2);
  y = zeros(1,numNeutrons_born*2);
  weight = ones(1,numNeutrons_born*2);
  iGroup = ones(1,numNeutrons_born*2);

%--------------------------------------------------------------------------
% Neutrons are assumed born randomly distributed in the cell with weight 1
% with sampled fission energy spectrum
  numNeutrons = numNeutrons_born; 
  for iNeutron = 1:numNeutrons
      x(iNeutron) = rand()*pitch;
      y(iNeutron) = rand()*pitch;
      weight(iNeutron) = 1;
    % Sample the neutron energy group
      iGroup(iNeutron) = find(cumsum(fuel.chi) >= rand(), 1, 'first');
  end

%--------------------------------------------------------------------------
% Prepare vectors for keff and standard deviation of keff
  keff_expected = ones(1,numCycles_active);
  sigma_keff = zeros(1,numCycles_active);
  keff_active_cycle = ones(1,numCycles_active);
  virtualCollision = false;

% Main (power) iteration loop
  for iCycle = 1:(numCycles_inactive + numCycles_active)

    % Normalize the weights of the neutrons to make the total weight equal to
    % numNeutrons_born (equivalent to division by keff_cycle)
      weight = (weight ./ sum(weight,2)) * numNeutrons_born;
      weight0 = weight;

    %----------------------------------------------------------------------
    % Loop over neutrons
      for iNeutron = 1:numNeutrons

          absorbed = false;

        %------------------------------------------------------------------
        % Neutron random walk cycle: from emission to absorption

          while ~absorbed

           % Sample free path length according to the Woodcock method
             freePath = -log(rand())/SigTmax(iGroup(iNeutron));

             if ~virtualCollision				
               % Sample the direction of neutron flight assuming both
               % fission and scattering are isotropic in the lab (a strong
               % assumption!)
                 teta = pi*rand();
                 phi = 2.0*pi*rand();
                 dirX = sin(teta)*cos(phi);
                 dirY = sin(teta)*sin(phi);
             end

           % Fly
             x(iNeutron) = x(iNeutron) + freePath * dirX;
             y(iNeutron) = y(iNeutron) + freePath * dirY;             

           % If outside the cell, find the corresponding point inside the
           % cell
             while x(iNeutron) < 0, x(iNeutron) = x(iNeutron) + pitch; end
             while y(iNeutron) < 0, y(iNeutron) = y(iNeutron) + pitch; end
             while x(iNeutron) > pitch, x(iNeutron) = x(iNeutron) - pitch; end
             while y(iNeutron) > pitch, y(iNeutron) = y(iNeutron) - pitch; end

           % Find the total and scattering cross sections
             if x(iNeutron) > 0.9 && x(iNeutron) < 2.7                     % INPUT
                SigA = fuel.SigF(iGroup(iNeutron)) + fuel.SigC(iGroup(iNeutron)) + fuel.SigL(iGroup(iNeutron));
                SigS = fuel.SigS{1+0}(iGroup(iNeutron),:)';
                SigP = fuel.SigP(iGroup(iNeutron));
             elseif x(iNeutron) < 0.7 || x(iNeutron) > 2.9                 % INPUT
                SigA = cool.SigC(iGroup(iNeutron)) + cool.SigL(iGroup(iNeutron));
                SigS = cool.SigS{1+0}(iGroup(iNeutron),:)';
                SigP = 0;
             else
                SigA = clad.SigC(iGroup(iNeutron)) + clad.SigL(iGroup(iNeutron));
                SigS = clad.SigS{1+0}(iGroup(iNeutron),:)';
                SigP = 0;
             end
             
           % Find the other cross sections ...
           % ... scattering
             SigS_sum = sum(SigS);
           % ... total
             SigT = SigA + SigS_sum;
           % ... virtual
             SigV = SigTmax(iGroup(iNeutron)) - SigT;
             
           % Sample the type of the collision: virtual (do nothing) or real
             if SigV/SigTmax(iGroup(iNeutron)) >= rand() % virtual collision

                virtualCollision = true;

             else % real collision

                virtualCollision = false;

              % Sample type of the collision: scattering or absorption
                if SigS_sum/SigT >= rand() % isotropic scattering

                 % Score scatterings with account for weight divided by the
                 % total scattering cross section
                   detectS(iGroup(iNeutron)) = detectS(iGroup(iNeutron)) + weight(iNeutron)/SigS_sum;

                 % Sample the energy group of the secondary neutron
                   iGroup(iNeutron) = find(cumsum(SigS)/SigS_sum >= rand(), 1, 'first');

                else % absorption

                   absorbed = true;

                 % Neutron is converted to the new fission neutron with
                 % the weight increased by eta
                   weight(iNeutron) = weight(iNeutron) * (SigP/SigA);
                   
                 % Sample the energy group for the new-born neutron
                   iGroup(iNeutron) = find(cumsum(fuel.chi) >= rand(), 1, 'first');

                end % scattering or absorption
             end % virtual or real
          end % of neutron random walk cycle: from emission to absorption
      end % of loop over neutrons

    %----------------------------------------------------------------------
    % Russian roulette
      for iNeutron = 1:numNeutrons
          terminateP = 1 - weight(iNeutron)/weight0(iNeutron);          
          if terminateP >= rand()
             weight(iNeutron) = 0; % killed
          elseif terminateP > 0
             weight(iNeutron) = weight0(iNeutron); % restore the weight
          end
      end

    %----------------------------------------------------------------------
    % Clean up absorbed or killed neutrons

      x(weight == 0) = [];
      y(weight == 0) = [];
      iGroup(weight == 0) = [];     
      weight(weight == 0) = [];
      numNeutrons = size(weight,2);

    %----------------------------------------------------------------------
    % Split too "heavy" neutrons

      numNew = 0;
      for iNeutron = 1:numNeutrons
          if weight(iNeutron) > 1
           % Truncated integer value of the neutron weight
             N = floor(weight(iNeutron));
           % Sample the number of split neutrons
             if weight(iNeutron)-N > rand(), N = N + 1; end
           % Change the weight of the split neutron
             weight(iNeutron) = weight(iNeutron)/N;
           % Introduce new neutrons
             for iNew = 1:N-1
                 numNew = numNew + 1;
                 x(numNeutrons + numNew) = x(iNeutron);
                 y(numNeutrons + numNew) = y(iNeutron);
                 weight(numNeutrons + numNew) = weight(iNeutron);
                 iGroup(numNeutrons + numNew) = iGroup(iNeutron);
             end
          end
      end
    % Increase the number of neutrons
      numNeutrons = numNeutrons + numNew;

    %----------------------------------------------------------------------
    % k-eff in a cycle equals the total weight of the new generation over
    % the total weight of the old generation (the old generation weight =
    % numNeutronsBorn)
      keff_cycle = sum(weight,2)/sum(weight0,2);
            
      iActive = iCycle - numCycles_inactive;
      if iActive <= 0
         fprintf('Inactive cycle = %3i/%3i; k-eff cycle = %8.5f; numNeutrons = %3i\n', ...
                 iCycle,numCycles_inactive,keff_cycle,numNeutrons);
      else
       % k-effective of the cycle
         keff_active_cycle(iActive) = keff_cycle;

       % k-effective of the problem
         keff_expected(iActive) = mean(keff_active_cycle(1:iActive));

       % Standard deviation of k-effective
         sigma_keff(iActive) = sqrt( sum( ( keff_active_cycle(1:iActive) - keff_expected(iActive) ).^2 ) ...
                                     / max(iActive-1,1) / iActive );

         fprintf('Active cycle = %3i/%3i; k-eff cycle = %8.5f; numNeutrons = %3i; k-eff expected = %9.5f; sigma = %9.5f\n', ...
                 iCycle-numCycles_inactive, numCycles_active, keff_cycle, numNeutrons, keff_expected(iActive), sigma_keff(iActive));
      end

  end % of main (power) iteration

%--------------------------------------------------------------------------
% Write down the results  
  fileID = fopen('resultsPWR.m','w');
     fprintf(fileID,'function s = resultsPWR\n\n');
     fprintf(fileID,'%% Results for 2D neutron transport calculation in the PWR-like unit cell using method of Monte Carlo\n\n');
     fprintf(fileID,'%% Number of source neutrons per k-eff cycle is %3i\n', numNeutrons_born);
     fprintf(fileID,'%% Number of inactive source cycles to skip before starting k-eff accumulation %3i\n', numCycles_inactive);
     fprintf(fileID,'%% Number of active source cycles for k-eff accumulation %3i\n\n', numCycles_active);

   % Stop stopwatch
     elapsedTime = toc;
     fprintf(fileID,'%% elapsedTime = %12.5e s / %12.5e min / %12.5e hrs;\n\n',elapsedTime, elapsedTime/60, elapsedTime/3600);
     fprintf(fileID,'s.keff_expected = %12.5f;\n',keff_expected(end));
     fprintf(fileID,'s.sigma = %12.5f;\n\n',sigma_keff(end));
     fprintf(fileID,'s.keffHistory = ['); fprintf(fileID,'%12.5e ',keff_expected); fprintf(fileID,'];\n');
     fprintf(fileID,'s.keffError = ['); fprintf(fileID,'%12.5e ',sigma_keff); fprintf(fileID,'];\n\n');

     fprintf(fileID,'%% Energy mesh, eV\n');
     fprintf(fileID,'s.eg = [');   fprintf(fileID,'%12.5e ',(fuel.eg(1:ng)+fuel.eg(2:ng+1))/2);   fprintf(fileID,'];\n\n');% group energies

     du = log( fuel.eg(2:ng+1)./fuel.eg(1:ng) )';
     flux_du = detectS' ./ du;
     fprintf(fileID,'%% Flux per unit lethargy, a.u.\n');
     fprintf(fileID,'s.flux = ['); fprintf(fileID,'%12.5e ',flux_du); fprintf(fileID,'];\n\n');


     fprintf(fileID,'end');
  fclose(fileID);

%--------------------------------------------------------------------------
% Plot the k-effective
  f = figure('visible','off');
  plot(keff_expected,'-r');
  hold on;
  plot(keff_expected+sigma_keff,'--b');
  hold on;
  plot(keff_expected-sigma_keff,'--b');
  grid on;
  xlabel('Iteration number');
  ylabel('k-effective');
  legend('k_{eff}','k_{eff} \pm \sigma')
  axis tight;
  saveas(f, 'MC_01_keff.pdf');

% Plot the spectrum  
  f = figure('visible','off');
  semilogx((fuel.eg(1:ng)+fuel.eg(2:ng+1))/2, flux_du);
  grid on;
  xlabel('Energy, eV');
  ylabel('Neutron flux per unit lethargy, a.u.');
  saveas(f, 'MC_02_flux_lethargy.pdf');

end % of function