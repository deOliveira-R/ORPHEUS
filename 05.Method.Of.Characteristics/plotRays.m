% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.

function plotRays(nNodes, delta, tit, fileName)

  MIN = 0;
  MAX = delta*(nNodes-1);
   
  f = figure('visible','off');
  axis equal;
  axis([MIN MAX -MAX MIN]);
  xlabel('Distance from the cell centre (cm)');
  ylabel('Distance from the cell centre (cm)');
  title(tit);


% horizontal lines
  X = [MIN MAX];
  Y = [MIN MIN];
  for i=1:nNodes
      line(X, Y, 'Color', 'k');
      Y = Y - delta;
  end
% vertical lines
  X = [MIN MIN];
  Y = [-MAX MIN];
  for i=1:nNodes
      line(X, Y, 'Color', 'k');
      X = X + delta;
  end
% 45-deg lines
  X = [MIN MIN+delta];
  Y = [MIN-delta MIN];
  for i=1:2*nNodes
      line(X, Y, 'Color', 'k');
      X(2) = X(2) + delta;
      Y(1) = Y(1) - delta;
  end
% 135-deg lines
  X = [MIN MIN+delta];
  Y = [-MAX+delta -MAX];
  for i=1:2*nNodes
      line(X, Y, 'Color', 'k');
      X(2) = X(2) + delta;
      Y(1) = Y(1) + delta;
  end
  saveas(f, fileName);

end