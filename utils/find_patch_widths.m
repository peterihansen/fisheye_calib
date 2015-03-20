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


function Pw = find_patch_widths(gpoints, i, imax, j, jmax, nr, nc)

    if i == 1
        i_index = [1];
    elseif i == imax
        i_index = [-1];
    else 
        i_index = [-1 1];
    end
    
    if j == 1
        j_index = [1];
    elseif j == jmax
        j_index = [-1];
    else 
        j_index = [-1 1];
    end

    c = 0;
    for a = 1:length(i_index)
        for b = 1:length(j_index)
           c = c + 1;
           side(c) = sqrt( ...
               (gpoints(i,j,1)-gpoints(i+i_index(a),j+j_index(b),1))^2 + ...
               (gpoints(i,j,2)-gpoints(i+i_index(a),j+j_index(b),2))^2 );
        end
    end
    pwidth = round(min(side)/6);    % 3
    
    Pw(1:4) = pwidth;
    x = round(gpoints(i,j,1));
    y = round(gpoints(i,j,2));
    if (x-pwidth) < 1
        Pw(1) = x -1;
    end
    if (x+pwidth) > nc
        Pw(2) = nc - x;
    end
    if (y-pwidth) < 1
        Pw(3) = y - 1;
    end
    if (y+pwidth) > nr
        Pw(4) = nr - y;
    end
    
    