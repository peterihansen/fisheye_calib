%
%   cam = set_unified_model(cam)
%
%   Finds an approximate estimate of the camera model
%
%   cam.pp = [u;v]  principal point estimate   (input)
%   cam.r90 = radius from pp for angle colatitude = 90deg.
%   
%   cam.l = the unified point of projection
%   cam.m = the distance of image plane from center of view sphere
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

function cam = set_unified_model(cam)

    % Point of projection (l = 2.0)
    % This should be in the ballpark for fisheye cameras
    % (0 = persp, 0->1 = hyperbolic cata, 1 = parabolic cata, >1 = most fisheye) 
    cam.l = 1.5;
    
    % Distance of image place from center of view sphere
    cam.m = cam.l * (cam.r90 - 1); 
end
    
    

