% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function is called by the iterative solver of the linear algebraic 
% equations and return Ax when x is an input

function ax = funCORE1D(solution)

  global nNodes dz Sig D_ nGroups J;

  fi = reshape(solution,nGroups,nNodes);

% Evaluate the gradient by a simple finite difference scheme
  dfidz = [(fi(:,1)-0)/(0.5*dz), diff(fi,1,2)/dz, (0-fi(:,nNodes))/(0.5*dz)];
  
% Fick's law 
  J = -D_ .* dfidz;
  
  Ax = diff(J,1,2)/dz + (Sig.A + Sig.S) .* fi - flipud(Sig.S) .* flipud(fi);
 
% Make one column-vector  
  ax = reshape(Ax,[],1);
end