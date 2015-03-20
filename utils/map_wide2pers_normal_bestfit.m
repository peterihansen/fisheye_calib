%====================================================================
% function im = map_wide2pers_normal_bestfit
%
% Maps a wide angle image to a perspective image using the unified
% image model.  The mapping is controlled by the focal length F, point
% of projection L, image centres cX and cY, scaling parameter 
% sX and sY, and a rotation of the sphere given by two euler angles.  
% The size of the output image is the same as the input.
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

function imvars = map_wide2pers_normal_bestfit(im, cam, imvars, ah)           

    title(ah(1),'Generating New Image','FontWeight','Bold','FontSize',12);                                                      
    drawnow
    
    imvars.remapped = imvars.remapped + 1;
                                            
    %----------------------------------------------------------
    % The data is the points on the sphere (circular indexing)    
    %----------------------------------------------------------     
    [nr,nc] = size(im);
    nr2 = 2*nr;
    nc2 = 2*nc;
    
    % Get Euler angles for the normal to the plane
    beta = imvars.normal(1);
    gamma = imvars.normal(2);
    alpha = imvars.normal(3);
    
    % Set up the rotation matrices (defined by normal to the plane)
    Ry = euler_matrix('y',beta);
    Rz = euler_matrix('z',gamma);
    Rz2 = euler_matrix('z',alpha);
    
    % Find the position on the perspective image plane of the corner
    % points of the border
    k = 0;
    l = length(imvars.border_wide{1}{1}(:,1));
    R = Rz * Ry * Rz2;
    for j = 1:2
        for i = 1:l-1:l
            k = k + 1;
            S = cam.map_img2ray(cam, imvars.border_wide{1}{j}(i,:)');
            S = R' * S;
            x_border(1,k) = S(1)/S(3);
            y_border(1,k) = S(2)/S(3);               
        end
    end

    
    % Use the min and max values to optimise the view of the grid
    % and set a suitable pinhole matrix K
    %-----------------------------------------------------------
    minX = min(x_border);   %round(min(x_border));
    maxX = max(x_border);   %round(max(x_border));
    minY = min(y_border);   %round(min(y_border));
    maxY = max(y_border);   %round(max(y_border));
    dX = maxX - minX;
    dY = maxY - minY;
    fX = nc2 / (1.2*dX);
    fY = nr2 / (1.2*dY);
    f = min(fX,fY);
    u0 = (nc2/2) - (f * (minX + maxX)/2);
    v0 = (nr2/2) - (f * (minY + maxY)/2);
    K = [f 0 u0; 0 f v0; 0 0 1];
    
    
    % Make the new orthonormal image
    %-----------------------------------------------------------
    [u,v] = meshgrid((1:nc2)-1, (1:nr2)-1);
    S = inv(K) * [u(:)'; v(:)'; ones(1,nr2*nc2)];
    S = R * (S ./ repmat(sqrt(sum(S.^2)),[3 1]));
    
    U_orig = cam.map_ray2img(cam, S);
    u = reshape(U_orig(1,:),[nr2 nc2])+1;
    v = reshape(U_orig(2,:),[nr2 nc2])+1;
    im_wide = interp2(im, u, v, 'cubic');

    % Set the outputs
    imvars.remap_img = conv2(im_wide, fspecial('gaussian',3,1), 'same');
    imvars.remap_K = K;
end

