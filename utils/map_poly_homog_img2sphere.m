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


%function [varargout] = map_poly_img2sphere(U, H, Cxy, pvars)
function [varargout] = map_poly_img2sphere(U, H, pvars)

    % Make homogeneous
    U = [U; 1];

    % Remove camera center and homography
    U = inv(H) * (U);% - [Cxy;0]);
    U = U(1:2) ./ repmat(U(3),[2 1]);
    
    % Polar coords
    phi = atan2(U(2),U(1));
    r = sqrt(U(1)^2 + U(2)^2);
     
    % Polynomial roots
    proots = roots([flipud(pvars); -r]);
    ind = find(proots >=0 & imag(proots) == 0);
    
    if isempty(ind)
        ind = find(imag(proots) ~= 0);
        theta = abs(proots(ind(1)));
    else
        if length(ind) == 1
            theta = proots(ind);
        else
            theta = min(proots(ind));
        end
    end
    
    sx = sin(theta) * cos(phi);
    sy = sin(theta) * sin(phi);
    sz = cos(theta);
        
    if nargout == 1
        varargout = {[sx;sy;sz]};
    else
        varargout(1) = {sx};
        varargout(2) = {sy};
        varargout(3) = {sz};
    end
end
    