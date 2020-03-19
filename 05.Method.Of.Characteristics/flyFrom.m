% ---------------------------------------------------------
% Matlab-based Open-source Reactor Physics Education System
% ---------------------------------------------------------
% Author: Konstantin Mikityuk, Paul Scherrer Institute, 2015-2017.
%
% Given coordinate ix, iy the function flyFrom updates angular flux in
% the neighbouring node according to the direction of the neutron fly

function flyFrom(ix, iy, direction)

  global g fi SigT q d

  switch direction
      
      case d.EAST
          if ix == g.nNodes
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.EAST), SigT{ix,iy}', g.delta/2, q{ix,iy});
          fi{ix+1,iy}(:,d.EAST) = fly(fi_, SigT{ix+1,iy}', g.delta/2, q{ix+1,iy});
          
      case d.NORTH_EAST
          if ix == g.nNodes || iy == 1
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.NORTH_EAST), SigT{ix,iy}', g.delta/sqrt(2), q{ix,iy});
          fi{ix+1,iy-1}(:,d.NORTH_EAST) = fly(fi_, SigT{ix+1,iy-1}', g.delta/sqrt(2), q{ix+1,iy-1});
          
      case d.NORTH
          if iy == 1
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.NORTH), SigT{ix,iy}', g.delta/2, q{ix,iy});
          fi{ix,iy-1}(:,d.NORTH) = fly(fi_, SigT{ix,iy-1}', g.delta/2, q{ix,iy-1});
          
      case d.NORTH_WEST
          if ix == 1 || iy == 1
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.NORTH_WEST), SigT{ix,iy}', g.delta/sqrt(2), q{ix,iy});
          fi{ix-1,iy-1}(:,d.NORTH_WEST) = fly(fi_, SigT{ix-1,iy-1}', g.delta/sqrt(2), q{ix-1,iy-1});
          
      case d.WEST
          if ix == 1
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.WEST), SigT{ix,iy}', g.delta/2, q{ix,iy});
          fi{ix-1,iy}(:,d.WEST) = fly(fi_, SigT{ix-1,iy}', g.delta/2, q{ix-1,iy});
          
      case d.SOUTH_WEST
          if ix == 1 || iy == g.nNodes
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.SOUTH_WEST), SigT{ix,iy}', g.delta/sqrt(2), q{ix,iy});
          fi{ix-1,iy+1}(:,d.SOUTH_WEST) = fly(fi_, SigT{ix-1,iy+1}', g.delta/sqrt(2), q{ix-1,iy+1});
          
      case d.SOUTH
          if iy == g.nNodes
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.SOUTH), SigT{ix,iy}', g.delta/2, q{ix,iy});
          fi{ix,iy+1}(:,d.SOUTH) = fly(fi_, SigT{ix,iy+1}', g.delta/2, q{ix,iy+1});
          
      case d.SOUTH_EAST
          if ix == g.nNodes || iy == g.nNodes
             return
          end
          fi_ = fly(fi{ix,iy}(:,d.SOUTH_EAST), SigT{ix,iy}', g.delta/sqrt(2), q{ix,iy});
          fi{ix+1,iy+1}(:,d.SOUTH_EAST) = fly(fi_, SigT{ix+1,iy+1}', g.delta/sqrt(2), q{ix+1,iy+1});
  end
end
