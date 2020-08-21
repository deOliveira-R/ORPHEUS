% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function homogeneousAqueousReactor calculates the neutron spectrum
% and k-infinitive by solving an eigenvalue problem for the homogeneous
% aqueous reactor.

function homogeneousAqueousReactor
  
% Number of energy groups:
  ng = 421;

% Path to macroscopic cross section data:
  path(path,['..' filesep '02.Macro.XS.421g']);

% Fill structure s with the data:
  s = macro421_H2OU__294K;

% Set initial values for flux 
  fi = ones(ng,1);

% Main iteration loop (5 iterations should be enough to converge)
  for i=1:5
   % pRate is total neutron production rate, a scalar
     pRate = (s.SigP + 2*sum(s.Sig2,2)') * fi;
   % aRate is total neutron absorption rate, a scalar
     aRate = (s.SigC + s.SigF + s.SigL + sum(s.Sig2,2)') * fi;
     
   % We evaluate the multiplication factor as a ratio of neutron
   % production rate and neutron absorption rate (there is no neutron
   % leakage in the infinite reactor):
     kInf = pRate / aRate;
     fprintf('kInf = %9.5f\n',kInf);   

   % Slowing down equation: diag(SigT).*fi = chi' * pRate / kInf +
   % SigS0'.*fi + 2*Sig2'.*fi
   % Flux:
     fi = inv(diag(s.SigT) - s.SigS{1+0}' - 2*s.Sig2') * (s.chi' * pRate / kInf);        %#ok<*MINV>
   
   % Normalize the flux so that the total production rate is 100 n/s per
   % cm3
     fi = fi / pRate * 100;
  end
  
% Write down the results  
  fileID = fopen('resultsAqueousReactor.m','w');
     fprintf(fileID,'function s = resultsAqueousReactor\n\n');
     fprintf(fileID,'%% Results for neutron slowing down in an infinite homogeneous aqueous reactor\n\n');
     
     fprintf(fileID,'%% k-Infinity\n');
     fprintf(fileID,'s.kInf = %9.5f;\n\n',kInf);
     
     fprintf(fileID,'%% Energy mesh, eV\n');
     fprintf(fileID,'s.eg = [');   fprintf(fileID,'%12.5e ',(s.eg(1:ng)+s.eg(2:ng+1))/2);   fprintf(fileID,'];\n\n');% group energies
     
     de = (s.eg(2:ng+1) - s.eg(1:ng))';
     fprintf(fileID,'%% Energy bins, eV\n');
     fprintf(fileID,'s.de = [');   fprintf(fileID,'%12.5e ',de);   fprintf(fileID,'];\n\n');
     
     fprintf(fileID,'%% Flux per unit energy, 1/(cm2*s*eV)\n');
     fprintf(fileID,'s.fi_de = [');  fprintf(fileID,'%12.5e ',fi(1:ng)./de);  fprintf(fileID,'];\n\n');
     
     du = log( s.eg(2:ng+1)./s.eg(1:ng) )';
     fprintf(fileID,'%% Lethargy bins\n');
     fprintf(fileID,'s.du = [');   fprintf(fileID,'%12.5e ',du);   fprintf(fileID,'];\n\n');
     
     fprintf(fileID,'%% Flux per unit lethargy, 1/(cm2*s)\n');
     fprintf(fileID,'s.fi_du = [');  fprintf(fileID,'%12.5e ',fi(1:ng) ./ du);  fprintf(fileID,'];\n\n');
     
     fprintf(fileID,'%% Total flux, 1/(cm2*s)\n');
     fprintf(fileID,'s.fi_tot = %12.5e;\n\n',sum(fi));
     
     fprintf(fileID,'%% One-group production cross section, 1/cm\n');
     fprintf(fileID,'s.sigP = %12.5e;\n\n',pRate/sum(fi));
     
     fprintf(fileID,'%% One-group absorption cross section, 1/cm\n');
     fprintf(fileID,'s.sigA = %12.5e;\n\n',aRate/sum(fi));
     
     fprintf(fileID,'end');
  fclose(fileID);
  
% Plot the spectrum  
  f = figure('visible','off');
  loglog((s.eg(1:ng)+s.eg(2:ng+1))/2, fi(1:ng)./de);
  grid on;
  xlabel('Energy, eV');
  ylabel('Neutron flux per unit energy, 1/(cm2*s*eV)');
  saveas(f, 'AQ_01_homogeneousAqueousReactor_flux_energy.pdf');

  f = figure('visible','off');
  semilogx((s.eg(1:ng)+s.eg(2:ng+1))/2, fi(1:ng) ./ du);
  grid on;
  xlabel('Energy, eV');
  ylabel('Neutron flux per unit lethargy, 1/(cm2*s)');
  saveas(f, 'AQ_02_homogeneousAqueousReactor_flux_lethargy.pdf');
  
end