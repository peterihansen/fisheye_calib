%====================================================================
% function P = map_sphere2sphere(S, normal, dir)
%
% Maps point x,y,z on a unit radius sphere to the image plane. The 
% mapping is  given by the focal length F and point of projection L 
% of the sphere
% 
% Returns the point on the image plane P
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

function [varargout] = map_sphere2sphere(S, normal, dir)

%     if dir == 1
%         beta = -normal(1);
%     else
%         beta = normal(1);
%     end
%     
%     gamma = normal(1);
%     Ry = euler_matrix('y',beta);
%     Rz = euler_matrix('z',gamma);
%     R = Rz * Ry;
%     
%     Sout = R * S';
%     
%     if nargout == 1
%         varargout = {Sout};
%     else
%         varargout(1) = {Sout(1)};
%         varargout(2) = {Sout(2)};
%         varargout(3) = {Sout(3)};
%     end

    
    beta = normal(1);
    gamma = normal(2);
    alpha = normal(3);
    Ry = euler_matrix('y',beta);
    Rz = euler_matrix('z',gamma);
    Rz2 = euler_matrix('z',alpha);
    R = Rz * Ry * Rz2;
    
    if dir == 0
        Sout = R * S;
    else
        Sout = R' * S;
    end
    
    
    if nargout == 1
        varargout = {Sout};
    else
        varargout(1) = {Sout(1)};
        varargout(2) = {Sout(2)};
        varargout(3) = {Sout(3)};
    end