%
%
%   This work was made possible by YSREP grant #1-019-2-008 from the Qatar 
%   National Research Fund (a member of Qatar Foundation).  
%
%   Copyright (C) 2015 Carnegie Mellon University (CMU)
%
%   This file is part of fisheye_calib
%
%   fisheye_calib is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%


function points = update_grid(gpoints, i, imax, j, jmax)

    % Set the new corner points of the grid 
    if (i == imax) && (j == jmax)
        points = gpoints;
        return
    end
    
    if i == imax
        corners(1,:) = [gpoints(i,j,1) gpoints(i,j,2)];
        corners(2,:) = [gpoints(i,jmax,1) gpoints(i,jmax,2)];
        plines{1}{1} = corners;
        num_grid = [jmax-j];
    elseif j == jmax
        corners(1,:) = [gpoints(i,j,1) gpoints(i,j,2)];
        corners(2,:) = [gpoints(imax,j,1) gpoints(imax,j,2)];
        plines{1}{1} = corners;
        num_grid = [imax-i];
    else
        corners(1,:) = [gpoints(i,j,1) gpoints(i,j,2)];
        corners(2,:) = [gpoints(imax,j,1) gpoints(imax,j,2)];
        corners(3,:) = [gpoints(imax,jmax,1) gpoints(imax,jmax,2)];
        corners(4,:) = [gpoints(i,jmax,1) gpoints(i,jmax,2)];
        plines = parallel_lines(corners);
        num_grid = [imax-i jmax-j]; % This must be in the right order
    end

    
    % Get the grid coords within this border
    minigrid_temp = grid_coords(plines, num_grid, 1); 
    
    
    if i == imax
        minigrid(:,:,1) = minigrid_temp(:,:,1)';
        minigrid(:,:,2) = minigrid_temp(:,:,2)';
    else
        minigrid = minigrid_temp;
    end

    
    points = gpoints;
    ai = 0;
    for a = i:imax
        ai = ai + 1;
        bj = 0;
        for b = j:jmax
            bj = bj + 1;
            points(a,b,1) = minigrid(ai,bj,1);
            points(a,b,2) = minigrid(ai,bj,2);
        end
    end
            