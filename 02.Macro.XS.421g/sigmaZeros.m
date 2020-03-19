% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% Given sigTtab -- the tables for every isotope of total cross sections 
% versus sigma-zeros for every energy group, the base points sig0tab, 
% the vector of atomic (number) densities, aDen, and the escape cross
% section, SigEscape, the function calculates sigma-zeros -- the background
% cross sections.
%
% Input:
%
% sigTtab -- a cell array containing nIso (number of isotopes in the
% mixture) 2D matrices of total microscopic cross sections of dimension
% nsigz x ng, where nsigz is the number of the base points sigma-zeros
% (could be different for different isotopes) and ng = 421 is the number of
% the energy groups. 
%
% sig0tab -- a cell array containg ng vectors of the base point of
% tabulated sigma-zeros.
%
% aDen -- a vector of length nIso with atomic densities of isotopes.
%
% SigEscape -- escape cross section (1/cm), for simple convex objects (such
% as plates, spheres, or cylinders) is given by S/(4V), where V and S are
% the volume and surface area of the object, respectively
%
% Output:
%
% sig0 -- 2D matrix of sigma-zeros of dimension (nIso x ng)
%
% sigT -- 2D matrix of total macroscopic cross sections corrected with 
% account for sigma-zero (nIso x ng)


function sig0 = sigmaZeros(sigTtab, sig0tab, aDen, SigEscape)

% Number of energy groups
  ng = 421;

% Define number of isotopes in the mixture
  nIso = length(aDen);

% first guess for sigma-zeros is 1e10 (infinite dilution)
  sig0 = 1e10 * ones(nIso,ng);

% Loop over energy group 
  for ig = 1:ng

    % Error to control sigma-zero iterations
      err = 1e10;
      nIter = 0;
 
    % sigma-sero iterations until the error is below selected tolerance (1e-6)
      while err > 1e-6

       % Loop over isotopes
         for iIso = 1:nIso
           % Find cross section for the current sigma-zero by interpolating
           % in the table
             if length(sig0tab{iIso}) == 1
                sigT(iIso,ig) = sigTtab{iIso}(1,ig);
             else
                log10sig0 = min(10,max(0,log10(sig0(iIso,ig))));
                sigT(iIso,ig) = interp1(log10(sig0tab{iIso})', sigTtab{iIso}(:,ig), log10sig0);
             end
         end

          err = 0;
        % Loop over isotopes
          for iIso = 1:nIso
            % Find the total macroscopic cross section for the mixture of
            % the background isotopes
              sum = 0;              
            % Loop over background isotopes
              for jIso = 1:nIso
                  if jIso ~= iIso
                     sum = sum + sigT(jIso,ig)*aDen(jIso);
                  end
              end
              tmp = (SigEscape + sum) / aDen(iIso);
              err = err + ( 1 - tmp/sig0(iIso,ig) )^2;
              sig0(iIso,ig) = tmp;
          end
          err = sqrt(err);
          nIter = nIter + 1;
          if nIter > 100
             fprintf('Error: too many sigma-zero iterations.\n')
             return
          end
      end
  end

end
