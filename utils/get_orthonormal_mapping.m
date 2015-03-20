%
%   Finds the normal to the calibration plane
%
%   This is major overkill for 4 points: this function was
%   originally used as part of parallel grid detection and calibration
%   with more than 4 points.
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

function imvars = get_orthonormal_mapping(imvars, cam, ah)
    
    title(ah(1),'Calibrating Normal to Plane',...
                'FontWeight','Bold', 'FontSize',14); 
    drawnow
    
    % Set the optimisation options 
    options = optimset('fminsearch');
    options.Dipslay = 'off';
    options.MaxFunEvals = 2000;
        
    % Map the points to the sphere
    for k = 1:length(imvars.border_wide)
        imvars.lineset = k;
        for i = 1:length(imvars.border_wide{k})
            imvars.line = i;
            S = cam.map_img2ray(cam, imvars.border_wide{k}{i}')';
            
            %u = imvars.border_wide{k}{i}(:,1) - cam.Cx;
            %v = imvars.border_wide{k}{i}(:,2) - cam.Cy;
            %clear S
            %for j = 1:length(u)
            %    S(j,:) = map_unified_img2sphere(u(j), v(j), cam.m, cam.l);
            %end
            imvars.border_sphere{k}{i} = S;
            
            % Fit great circles
            estimate_border_euler = imvars.border_euler{k}{i};
            [border_euler{k}{i},fval,exitflag] = fminsearch ...
                            (@objfun_euler_2param, estimate_border_euler);
        end
        % Get the points of intersection
        intercepts(k,:) = gcircle_intercepts(border_euler{k});
        intercepts(k+2,:) = -intercepts(k,:);
    end
    imvars.border_intercepts = intercepts;
        
    % Fit a great circle to the intercepts (antipodal points)
    estimate_border_normal = imvars.normal;
    [border_normal,fval,exitflag] = fminsearch ...
                        (@objfun_border_normal, estimate_border_normal);
    
    imvars.border_euler = border_euler;
    imvars.normal = border_normal;
        
    while border_normal(1) > pi
        border_normal(1) = border_normal(1) - 2*pi;
    end
    while border_normal(1) < -pi
        border_normal(1) = border_normal(1) + 2*pi;
    end
    imvars.normal = border_normal;
    
    % Find the final normal rotation
    points = imvars.border_intercepts;
    Ry = euler_matrix('y',border_normal(1));
    Rz = euler_matrix('z',border_normal(2));
    R = (Rz*Ry)';
    S1 = R * points(1,:)';
    %S2 = R * points(2,:)';
    imvars.normal(3) = atan2(S1(2), S1(1));
    
    
    % Objective function - Euler 2 param
    function f = objfun_euler_2param(angles)
    
        % The objective function is based on the radial distance
        % between the data points and the perpendicular point
        % on a circle of best fit.  The variables are the euler angles.
        points = imvars.border_sphere{imvars.lineset}{imvars.line};
   
        for i = 1:length(points(:,1))    
            rad_dist(i) = asin ...
                    (points(i,1) * sin(angles(1)) * cos(angles(2)) + ...
                      points(i,2) * sin(angles(1)) * sin(angles(2)) + ...
                      points(i,3) * cos(angles(1)));
        end
        f = sum(rad_dist.^2);
    end
        
    
    % Objective function - Border normal
    function f = objfun_border_normal(angles)
    
        % The objective function is based on the radial distance
        % between the data points and the perpendicular point
        % on a circle of best fit.  The variables are the euler angles.
    
        points = imvars.border_intercepts;
            
        for i = 1:length(points(:,1))    
            rad_dist(i) = asin ...
                   (points(i,1) * sin(angles(1)) * cos(angles(2)) + ...
                      points(i,2) * sin(angles(1)) * sin(angles(2)) + ...
                      points(i,3) * cos(angles(1)));
        end
        f = sum(rad_dist.^2);
    end
end     
             