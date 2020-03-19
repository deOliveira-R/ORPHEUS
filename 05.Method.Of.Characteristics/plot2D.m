% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function plots a 2D chart of function fun using the mesh for
% the 2D unit cell specified by delta and nNodes. The figure with the plot
% is saved as the fileName

function plot2D(nNodes, delta, fun, tit, fileName)

  f = figure('visible','off');
  axis equal;
  axis([0 delta*(nNodes-1) -delta*(nNodes-1) 0]);
  xlabel('Distance from the cell centre (cm)');
  ylabel('Distance from the cell centre (cm)');
  title(tit);

  fMax = 0;
  fMin = 1e20;
  for ix = 1:nNodes
      for iy = 1:nNodes
          fMax = max(fMax, fun(ix,iy));
          fMin = min(fMin, fun(ix,iy));
      end
  end

  y = 0;
  for iy = 1:nNodes
      x = 0;
      height = delta - delta/2 * (iy == 1 || iy == nNodes);
      y = y - height;
      for ix = 1:nNodes
          width = delta - delta/2 * (ix == 1 || ix == nNodes);
          color = (fun(ix,iy) - fMin)/(fMax - fMin);
          rectangle('Position',[x y width height],'FaceColor',[color 0 1-color],'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.2);
          x = x + width;
      end
  end

  saveas(f, fileName);
end
