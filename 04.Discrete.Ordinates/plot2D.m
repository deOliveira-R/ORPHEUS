% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% The function plots a 2D chart of function fun using the mesh specified 
% by nNodesX, nNodesY and delta. The figure with the plot is saved as the 
% fileName

function plot2D(nNodesX, nNodesY, delta, fun, tit, fileName)

  f = figure('visible','off');
  axis equal;
  axis([0 delta*(nNodesX-1) -delta*(nNodesY-1) 0]);
  xlabel('Distance from the cell centre (cm)');
  ylabel('Distance from the cell centre (cm)');
  title(tit);

  fMax = 0;
  fMin = 1e20;
  for ix = 1:nNodesX
      for iy = 1:nNodesY
          fMax = max(fMax, fun(ix,iy));
          fMin = min(fMin, fun(ix,iy));
      end
  end

  y = 0;
  for iy = 1:nNodesY
      x = 0;
      height = delta - delta/2 * (iy == 1 || iy == nNodesY) ;
      y = y - height;
      for ix = 1:nNodesX
          width = delta - delta/2 * (ix == 1 || ix == nNodesX) ;
          color = (fun(ix,iy) - fMin)/max(fMax - fMin);
          isnan(color);
          if isnan(color)
             color = 0;
          end
          rectangle('Position',[x y width height],'FaceColor',[color 0 1-color],'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.2);
          x = x + width;
      end
  end

  saveas(f, fileName);
end
