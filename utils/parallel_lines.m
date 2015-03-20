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

function plines = parallel_lines(data)

    % Set the parallel lines
    if length(data(:,1)) == 4
        lset{1} = [1 2; 4 3];
        lset{2} = [1 4; 2 3];
    else
        lset{1} = [1 2 3; 7 6 5];
        lset{2} = [1 8 7; 3 4 5];
    end
        
    % Get the output data in the required format
    for k = 1:length(lset)
        for j = 1:length(lset{k}(:,1))
            for i = 1:length(lset{k}(1,:))
                plines{k}{j}(i,:) = data(lset{k}(j,i),:);
            end
        end
    end