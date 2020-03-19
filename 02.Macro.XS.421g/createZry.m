% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function reads the MICROsopic group cross sections in the Matlab
% format and calculates from them the MACROscopic cross sections for
% natural mixture of zirconium isotopes which is the basis of zircalloy --
% fuel cladding material of the pressurized water reactor.

function createZry

% Global structure with geometry parameters
  global g

% number of energy groups
  Zry.ng = 421;
  
% Path to library:
  path(path,['..' filesep '00.Lib']);
  
% input and initialize the geometry of the PWR unit cell (the function is
% in '..\00.Lib')
  input_and_initialize_PWR_like;

% Zircaloy is composed of pure Zirconium (~1% of tin neglected). There are
% four stable izotopes of zirconium: Zr090, Zr091, Zr092, Zr094 and one
% very long-lived: Zr096 with the following molar fractions:
  molFrZr = [0.5145 0.1122 0.1715 0.1738 0.0280];

% Path to microscopic cross section data:
  path(path,['..' filesep '01.Micro.XS.421g']);

% Call the functions for Zr isotopes and store the data in the structures.
% As an example it is done below for 600K, change when other temperatures
% needed:
  Zr90 = micro_ZR090__600K;                                                % INPUT
  Zr91 = micro_ZR091__600K;                                                % INPUT
  Zr92 = micro_ZR092__600K;                                                % INPUT
  Zr94 = micro_ZR094__600K;                                                % INPUT
  Zr96 = micro_ZR096__600K;                                                % INPUT
  
  Zry.temp = 600;                                                          % INPUT
  Zry.eg = Zr90.eg; 
  
% Mass of one "average" Zr atom in atomic unit mass [a.u.m.]
  Zry.aw = Zr90.aw*molFrZr(1) + Zr91.aw*molFrZr(2) + Zr92.aw*molFrZr(3) + Zr94.aw*molFrZr(4) + Zr96.aw*molFrZr(5); 

% Path to material properties:
  path(path,['..' filesep '00.LIB']);
% The function returns material properties of Zry in structure clad
  [~, ~, clad] = matpro({}, {});
  
% Zircaloy density:
  Zry.den = clad.rho*1e-3;   % [g/cm3]
  rho = Zry.den*1.0e-24;     % [g/(barn*cm)]
  rho = rho / 1.660538e-24;  % [(a.u.m.)/(barn*cm)]
  rho = rho / Zry.aw;        % [number of Zr atoms/(barn*cm)]

% The names of isotopes
  Zry.isoName{1} = 'Zr90';
  Zry.isoName{2} = 'Zr91';
  Zry.isoName{3} = 'Zr92';
  Zry.isoName{4} = 'Zr94';
  Zry.isoName{5} = 'Zr96';

% The number densities of isotopes:
  Zry.numDen(1) = molFrZr(1) * rho;
  Zry.numDen(2) = molFrZr(2) * rho;
  Zry.numDen(3) = molFrZr(3) * rho;
  Zry.numDen(4) = molFrZr(4) * rho;
  Zry.numDen(5) = molFrZr(5) * rho;

% Prepare for sigma-zero iterations:
  sigTtab = {Zr90.sigT, Zr91.sigT, Zr92.sigT, Zr94.sigT, Zr96.sigT};
  sig0tab = {Zr90.sig0, Zr91.sig0, Zr92.sig0, Zr94.sig0, Zr96.sig0};
  aDen = Zry.numDen';

% SigEscape -- escape cross section, for simple convex objects (such as
% plates, spheres, or cylinders) is given by S/(4V), where V and S are the
% volume and surface area of the object, respectively
  SigEscape = 1 / (2*g.clad.rOut*100);
  
  fprintf('Sigma-zero iterations. ');
  Zry.sig0 = sigmaZeros(sigTtab, sig0tab, aDen, SigEscape);
  fprintf('Done.\n');

  fprintf('Interpolation of microscopic cross sections for the found sigma-zeros. ');
  sigCtab = {Zr90.sigC, Zr91.sigC, Zr92.sigC, Zr94.sigC, Zr96.sigC};
  sigLtab = {Zr90.sigL, Zr91.sigL, Zr92.sigL, Zr94.sigL, Zr96.sigL};
  for ig = 1:Zry.ng
    % Loop over isotopes
      for iIso = 1:5
        % Find cross sections for the sigma-zero
          if length(sig0tab{iIso}) == 1
             sigC(iIso,ig) = sigCtab{iIso}(1,ig);
             sigL(iIso,ig) = sigLtab{iIso}(1,ig);
          else
             log10sig0 = min(10,max(0,log10(Zry.sig0(iIso,ig))));
             sigC(iIso,ig) = interp1(log10(sig0tab{iIso})', sigCtab{iIso}(:,ig), log10sig0);
             sigL(iIso,ig) = interp1(log10(sig0tab{iIso})', sigLtab{iIso}(:,ig), log10sig0);
          end
      end
  end
  
  for j = 0:2
      sigS{1+j,1} = interpSigS(1+j, Zr90, Zry.sig0(1,:));
      sigS{1+j,2} = interpSigS(1+j, Zr91, Zry.sig0(2,:));
      sigS{1+j,3} = interpSigS(1+j, Zr92, Zry.sig0(3,:));
      sigS{1+j,4} = interpSigS(1+j, Zr94, Zry.sig0(4,:));
      sigS{1+j,5} = interpSigS(1+j, Zr96, Zry.sig0(5,:));
  end
  fprintf('Done.\n');
 
% Macroscopic cross section [1/cm] is microscopic cross section for the 
% "average" atom [barn] times the number density [number of atoms/(barn*cm)]
  Zry.SigC = sigC'*aDen;
  Zry.SigL = sigL'*aDen;
  for j = 0:2
      Zry.SigS{1+j} = sigS{1+j,1}*aDen(1) + sigS{1+j,2}*aDen(2) + sigS{1+j,3}*aDen(3) + sigS{1+j,4}*aDen(4) + sigS{1+j,5}*aDen(5);
  end
  Zry.Sig2 = Zr90.sig2*aDen(1) + Zr91.sig2*aDen(2) + Zr92.sig2*aDen(3) + Zr94.sig2*aDen(4) + Zr96.sig2*aDen(5);
  Zry.SigT = Zry.SigC + Zry.SigL + sum(Zry.SigS{1+0},2) + sum(Zry.Sig2,2);
  
% Number of fissile isotopes
  Zry.nFis = 0;
  
  Zry.SigP(1) = 0.0;

% Make a file name which include the isotope name and the temperature
  if Zry.temp < 1000
     matName = ['macro421_Zry__',num2str(round(Zry.temp)),'K']; % name of the file with a temperature index
  else
     matName = ['macro421_Zry_',num2str(round(Zry.temp)),'K']; % name of the file with a temperature index
  end
  
% Make a header for the file to be created with important parameters for
% which the macroscopic cross sections were generated
  Zry.header{1} = sprintf('%% ---------------------------------------------------------');
  Zry.header{2} = sprintf('%% Matlab-based Open-source Reactor Physics Education System');
  Zry.header{3} = sprintf('%% ---------------------------------------------------------');
  Zry.header{4} = sprintf('%% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.');
  Zry.header{5} = sprintf('%%');
  Zry.header{6} = sprintf('%% Macroscopic cross sections for zircaloy (assumed natural zirconium)');
  Zry.header{7} = sprintf('%% Zry temperature:    %6.1f K',Zry.temp);
  Zry.header{8} = sprintf('%% Zry density:        %6.3f g/cm3',Zry.den);
  
% Change the units of number density from 1/(barn*cm) to 1/cm2
  Zry.numDen = Zry.numDen*1e24;

% Finally create the file with macroscopic cross sections
  writeMacroXS(Zry,matName);
end