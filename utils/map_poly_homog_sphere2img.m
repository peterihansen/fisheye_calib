%
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

%function [varargout] = map_poly_sphere2img(S, H, Cxy, pvars)
function [varargout] = map_poly_sphere2img(S, H, pvars)

    % Get angular coordinates
    phi = atan2(S(2),S(1));
    theta = acos(S(3));
    
    % Radius in image
    r = 0;
    for i = 1:length(pvars)
        r = r + pvars(i)*theta^i;
    end

    
    % Add the homography and image center
    U = H * [r * cos(phi); r * sin(phi); 1];
    U = U(1:2) ./ [U(3); U(3)];
    %U = U + Cxy;
        
    if nargout == 1
        varargout = {U'};
    else
        varargout(1) = {U(1)};
        varargout(2) = {U(2)};
    end
end
    