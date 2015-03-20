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

function R = quat2rot(Q)
%R = quat2rot(q)
% R is a 3x3 rotation matrix from quaternion q
Q = permute(Q,[1 3 2]);

R = [
   Q(1,1,:).^2+Q(2,1,:).^2-Q(3,1,:).^2-Q(4,1,:).^2 2.0.*(Q(2,1,:).*Q(3,1,:)-Q(1,1,:).*Q(4,1,:))   2.0.*(Q(2,1,:).*Q(4,1,:)+Q(1,1,:).*Q(3,1,:))
   2.0.*(Q(2,1,:).*Q(3,1,:)+Q(1,1,:).*Q(4,1,:))   Q(1,1,:).^2-Q(2,1,:).^2+Q(3,1,:).^2-Q(4,1,:).^2 2.0.*(Q(3,1,:).*Q(4,1,:)-Q(1,1,:).*Q(2,1,:))
   2.0.*(Q(2,1,:).*Q(4,1,:)-Q(1,1,:).*Q(3,1,:))   2.0.*(Q(3,1,:).*Q(4,1,:)+Q(1,1,:).*Q(2,1,:))   Q(1,1,:).^2-Q(2,1,:).^2-Q(3,1,:).^2+Q(4,1,:).^2];

