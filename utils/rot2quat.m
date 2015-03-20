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


function q = rot2quat(R)
%q = rot2quat(R)

q0 = sqrt(R(1,1,:)+R(2,2,:)+R(3,3,:)+1)/2;

q1 = (R(3,2,:)-R(2,3,:))/4./q0;
q2 = (R(1,3,:)-R(3,1,:))/4./q0;
q3 = (R(2,1,:)-R(1,2,:))/4./q0;

i = find(q0 == 0);
if ~isempty(i)
   q1(i) = sqrt(abs(.5 *(R(2,2,i)+R(3,3,i))));
   q2(i) = sqrt(abs(.5 *(R(1,1,i)+R(3,3,i))));
   q3(i) = sqrt(abs(.5 *(R(1,1,i)+R(2,2,i))));
   
   j = i(find(q1(i) ~= 0));
   if ~isempty(j)
      q2(j) = q2(j) .* sign(R(1,2,j));
      q3(j) = q3(j) .* sign(R(1,3,j));
   end
   
   j = i(find(q1(i) == 0));
   if ~isempty(j)
      q3(j) = q3(j) .* sign(R(2,3,j));
   end
end
q = [q0(:)'; q1(:)'; q2(:)'; q3(:)'];