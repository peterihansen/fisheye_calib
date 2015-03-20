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


function intercepts = gcircle_intercepts(eulers)

    nconics = length(eulers);
    k = 0;
    
    % Get the Hessian normal for the plane defined by each great circle
    for i = 1:nconics-1
        
        % Get euler matrices
        %A1 = euler_matrix('x',eulers{i}(2));
        %B1 = euler_matrix('z',eulers{i}(1));
        
        Ry1 = euler_matrix('y',eulers{i}(1));
        Rz1 = euler_matrix('z',eulers{i}(2));
        
        % Get three points on the plane of the great circle
        x1 = [0; 0; 0];
        x2 = Rz1 * (Ry1 * [1;0;0]);
        x3 = Rz1 * (Ry1 * [0;1;0]);
        
        % Find the normal to the plane
        n_hat1 = cross(x3-x1,x2-x1);
        
        for j = i+1:nconics
            k = k + 1;
            
            % Get euler matrices             
            %A2 = euler_matrix('x',eulers{j}(2));
            %B2 = euler_matrix('z',eulers{j}(1));
            Ry2 = euler_matrix('y',eulers{j}(1));
            Rz2 = euler_matrix('z',eulers{j}(2));
            
            
            % Get three points on the plane of the great circle
            x1 = [0; 0; 0];  
            x2 = Rz2 * (Ry2 * [1;0;0]);
            x3 = Rz2 * (Ry2 * [0;1;0]);
        
            % Find the normal to the plane
            n_hat2 = cross(x3-x1,x2-x1);
          
            % Get line of intersection of the planes
            line_int = cross(n_hat1,n_hat2);
            
            % Find one of the points of intersection (+y value)
            if norm(line_int) ~= 0
                p1 = line_int / norm(line_int);
                p2 = -p1;
            
                if k == 1
                    p_ref = p1;
                    p_int = p1;
                else
                    if sum(abs(p1-p_ref)) < sqrt(2)
                        p_int = p1;
                    else
                        p_int = p2;
                    end
                end
            
                if sum(isnan(p_int)) == 0
                    intercepts(k,:) = p_int';
                else
                    k = k - 1;
                end
            else
                k = k - 1;
            end
        end
    end    
  