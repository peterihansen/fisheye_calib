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


function gpoints = grid_coords(plines, num_grid, divs_known)
    
    % If there is only one line, return the values along this line
    if length(plines) == 1
        diffxy = plines{1}{1}(2,:) - plines{1}{1}(1,:);
        m = diffxy(1,2) / diffxy(1,1);
        line_eq(1) = m;
        line_eq(2) = plines{1}{1}(1,2) - m*plines{1}{1}(1,1);
        Xspacing = (plines{1}{1}(2,1)-plines{1}{1}(1,1)) / num_grid;
        gpoints(:,:,1) = [plines{1}{1}(1,1):Xspacing:plines{1}{1}(2,1)]';
        gpoints(:,:,2) = gpoints(:,:,1) * line_eq(1) + line_eq(2); 
    else
                 
        % Determine the length of each line, and get the linear 
        % equation to the line (y=mx+c)
        for k = 1:2
            for i = 1:2
                a = length(plines{k}{i}(:,1));
                diffxy(i,:) = plines{k}{i}(a,:) - plines{k}{i}(1,:);
                m = diffxy(i,2) / diffxy(i,1);
                border_eq{k}{i}(1) = m;
                border_eq{k}{i}(2) = plines{k}{i}(1,2) - m*plines{k}{i}(1,1);
            end
            diffxy = diffxy.^2;
            side(k) = sum ( sqrt(diffxy(1,1) + diffxy(1,2)) + ...
                                        sqrt(diffxy(2,1) + diffxy(2,2)) );
        end

        
        % If the divisions along each line is unknown, work it out
        % based on the length of the lines
        if divs_known == 1;
            index = [1 2];
        else
            num_grid = sort(num_grid);
            [val,index] = sort(side);
        end
        
       
        % Find the points along each border line
        for k = 1:2
            for i = 1:2
                a = length(plines{k}{i}(:,1));
                
                % Find the spacing between the points
                Xspacing = (plines{k}{i}(a,1)-plines{k}{i}(1,1)) / ...
                                                     num_grid(index(k));
                % Get the points along the border line
                border_points{k}{i}(:,1) = ...
                    [plines{k}{i}(1,1):Xspacing:plines{k}{i}(a,1)];
                
                border_points{k}{i}(:,2) = border_points{k}{i}(:,1) * ...
                              border_eq{k}{i}(1) + border_eq{k}{i}(2);  
            end
        end
       
        
        % Get the equations to all the lines of the calibration grid
        for k = 1:2
            lines_eq{k}(:,1) = [ (border_points{k}{2}(:,2) - ...
                                  border_points{k}{1}(:,2)) ./ ...
                                 (border_points{k}{2}(:,1) - ...
                                  border_points{k}{1}(:,1)) ]; 
                          
            lines_eq{k}(:,2) = [ border_points{k}{1}(:,2) - ...
                             lines_eq{k}(:,1).*border_points{k}{1}(:,1) ];                         
        end
    
        gpoints = grid_line_intercepts(lines_eq);
    end