% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function returns the scattering matrix sigS which is a result of
% interpolation between scattering matrices for sigma-zero base points
% (s.SigS) for the target points Sig0 -- vector of length ng

function sigS = interpSigS(jLgn, s, Sig0)

% number of energy groups
  ng = 421;

% number of sigma-zeros
  nSig0 = length(s.sig0);
  if nSig0 == 1
     sigS = s.sigS{jLgn,1};
  else
     for iSig0 = 1:nSig0
         [ifrom, ito, tmp1(iSig0,:)] = find(s.sigS{jLgn,iSig0});
     end
     
%    number of non-zeros in a scattering matrix
     nNonZeros = size(tmp1,2);
     for i = 1:nNonZeros
         log10sig0 = min(10,max(0,log10(Sig0(ifrom(i)))));
         tmp2(i) = interp1(log10(s.sig0)', tmp1(:,i), log10sig0);
     end
     
     sigS = sparse(ifrom, ito, tmp2, ng, ng);
  end
  
end
