%
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

function points = grid_line_intercepts(lines_eq)

    for i = 1:length(lines_eq{1}(:,1))
        
        m1 = lines_eq{1}(i,1);
        c1 = lines_eq{1}(i,2);
        
        for j = 1:length(lines_eq{2}(:,1))
            
            m2 = lines_eq{2}(j,1);
            c2 = lines_eq{2}(j,2);
            
            points(i,j,1) = (c1 - c2) / (m2 - m1);
            points(i,j,2) = m1*points(i,j,1) + c1;
        end
    end