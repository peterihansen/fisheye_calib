% Euler matricies (USES RIGHT HAND RULE !!!)
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

function A = euler_matrix(axis,theta)

switch axis
    case 'x'
        A = [1      0           0;
             0      cos(theta)  -sin(theta);
             0      sin(theta)  cos(theta)];
         
    case 'y'
        A = [cos(theta)     0   sin(theta);
             0              1   0;
             -sin(theta)    0   cos(theta)];
         
    case 'z'
        A = [cos(theta)     -sin(theta) 0;
             sin(theta)     cos(theta)  0
             0              0           1];
         
    otherwise
        error('Enter an axis')
end