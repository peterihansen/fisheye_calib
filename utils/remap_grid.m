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

function points = remap_grid(gpoints, cam, imvars)

    [nr,nc,nd] = size(gpoints);
    u = gpoints(:,:,1);
    v = gpoints(:,:,2);
    S = inv(imvars.remap_K) * [u(:)'; v(:)'; ones(1,nr*nc)];
    S = S ./ repmat(sqrt(sum(S.^2)),[3 1]);
    S = map_sphere2sphere(S, imvars.normal, 0);
    U = cam.map_ray2img(cam, S);
    points(:,:,1) = reshape(U(1,:),[nr,nc]);
    points(:,:,2) = reshape(U(2,:),[nr,nc]);
end

%     for i = 1:length(gpoints(:,1,1))
%         for j = 1:length(gpoints(1,:,1))
%             
%             U = gpoints(i,j,1);
%             V = gpoints(i,j,2);
%             
%             % Correct for any offset of origin or rescaling
%             X = (U-1) * imvars.remap_spacing + imvars.remap_originX;
%             Y = (V-1) * imvars.remap_spacing + imvars.remap_originY;
%             
%             % Map point to sphere
%             Sn = map_unified_img2sphere(X, Y, imvars.m, 0);
%               
%             % De-correct normal
%             S(i,:) = map_sphere2sphere(Sn, imvars.normal, 0);
%                 
%             % Map to wide angle image
%             [x y] = map_unified_sphere2img(S(i,1), S(i,2), S(i,3), ...
%                                                    imvars.m, imvars.l);
%             
%             % Convert to image coordinates
%             points(i,j,1) = x + imvars.Cx;
%             points(i,j,2) = y + imvars.Cy;
%         end
%     end
    