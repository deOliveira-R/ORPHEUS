% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
%
% Convert 1D solution vector x to the array of angular fluxes fi

function fi = convert(solution)

  global g ng
  
% Define the cell array for angular fluxes
  fi = zeros(ng,g.N,g.nNodesX,g.nNodesY);
  
  flux = reshape(solution,ng,[]);
  nEq = 0;
  for iy = 1:g.nNodesY
      for ix = 1:g.nNodesX
          for n = 1:g.N
              if g.muZ(n) >= 0 && ~(ix==1 && g.muX(n)>0) && ~(ix==g.nNodesX && g.muX(n)<0) ...
			                   && ~(iy==1 && g.muY(n)>0) && ~(iy==g.nNodesY && g.muY(n)<0)
                 nEq = nEq + 1;
               % Unknowns are angular fluxes:
                 fi(:,n,ix,iy) = flux(:,nEq);
              end
          end
      end
  end
  for n = 1:g.N
      if g.muZ(n) < 0
         fi(:,n,:,:) = fi(:,g.nRefZ(n),:,:);
      end
  end

% Boundary conditions
  for n = 1:g.N
      if g.muX(n)>0; fi(:,n,1,:) = fi(:,g.nRefX(n),1,:); end;
      if g.muX(n)<0; fi(:,n,g.nNodesX,:) = fi(:,g.nRefX(n),g.nNodesX,:); end
  end
  for n = 1:g.N
      if g.muY(n)>0; fi(:,n,:,1) = fi(:,g.nRefY(n),:,1); end;
      if g.muY(n)<0; fi(:,n,:,g.nNodesY) = fi(:,g.nRefY(n),:,g.nNodesY); end
  end
end
