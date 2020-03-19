% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function is called by the iterative solver of the linear algebraic 
% equations and return Ax when x is an input

function ax = funDO(solution)

  global g SigT
  
% Convert 1D solution vector x to the cell array of angular flux fi
  fi = convert(solution);

  nEq = 0;
  for iy = 1:g.nNodesY
      for ix = 1:g.nNodesX
          for n = 1:g.N
              if g.muZ(n) >= 0 && ~(ix==1 && g.muX(n)>0) && ~(ix==g.nNodesX && g.muX(n)<0) ...
			                   && ~(iy==1 && g.muY(n)>0) && ~(iy==g.nNodesY && g.muY(n)<0)
               % Gradients
                 [dfidx, dfidy] = gradients(n,ix,iy,fi);
				 
                 nEq = nEq + 1;
                 LHS(:,nEq) = g.muX(n)*dfidx + g.muY(n)*dfidy +  SigT{ix,iy}' .* fi(:,n,ix,iy);
              end
          end % end of cycle over directions
      end % of cycle over y-axis
  end % of cycle over x-axis

% Make 1D vector
  ax = reshape(LHS,[],1);

end

% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2016.
%
% The function returns dfidx and dfidy which are gradients of the angular
% fluxes calculated according to the diamond scheme

function [dfidx, dfidy] = gradients(n,ix,iy,fi)

  global g

  if g.muX(n) > 0
     if ix == 1
        dfiX = fi(:,g.nRefX(n),ix,iy) - fi(:,g.nRefX(n),ix+1,iy);
     else %if ix > 1
        dfiX = fi(:,n,ix,iy) - fi(:,n,ix-1,iy);
     end
  else %if g.muX(n) <= 0
     if ix == g.nNodesX
        dfiX = fi(:,g.nRefX(n),ix-1,iy) - fi(:,g.nRefX(n),ix,iy);
     else % if ix < g.nNodesX
        dfiX = fi(:,n,ix+1,iy) - fi(:,n,ix,iy);
     end
  end % of g.muX(n)

  if g.muY(n) > 0
     if iy == 1
        dfiY = fi(:,g.nRefY(n),ix,iy) - fi(:,g.nRefY(n),ix,iy+1);
     else % if iy > 1
        dfiY = fi(:,n,ix,iy) - fi(:,n,ix,iy-1);
     end
  else % if g.muY(n) <= 0
     if iy == g.nNodesY
        dfiY = fi(:,g.nRefY(n),ix,iy-1) - fi(:,g.nRefY(n),ix,iy);
     else % if iy < g.nNodesY
        dfiY = fi(:,n,ix,iy+1) - fi(:,n,ix,iy);
     end
  end % of g.muY(n)

  dfidx = dfiX / g.delta;
  dfidy = dfiY / g.delta;
end