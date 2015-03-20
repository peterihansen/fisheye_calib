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

function [gpoints_original, Rnormals] = fisheye_grid_detect(img, cam, numgrid)

    % Initialise the figure
    ah = init_grid_detect_figure(size(img,1), size(img,2));
    
    % imgvars struct: needed for the remapping etc.
    mapvars.remapped = false;
    mapvars.remap_img = img;
    mapvars.remap_originX = 0;
    mapvars.remap_originY = 0;
    mapvars.remap_spacing = 0;
    mapvars.numgrid = numgrid;
    %mapvars.Cx = cam.Cx;
    %mapvars.Cy = cam.Cy;
    %mapvars.m = cam.m;
    %mapvars.l = cam.l;
    for k = 1:2
        for i = 1:2
            mapvars.border_wide{k}{i} = [0 0; 0 0];      % border corners
            mapvars.border_sphere{k}{i} = zeros(4,3);    % on sphere
            mapvars.border_euler{k}{i} = [0 0];          % euler angles
        end
    end
    mapvars.border_intercepts = [];  % Intercepts of the border lines
    mapvars.normal = [0; 0];          % the normal to the plane
    
    % Ugly, but quick hack to avoid all the loops and conditions
    % needed when using the 8 point selection etc.
    mapvars = fisheye_clickborder(cam, mapvars, ah, false);
    
    mapvars = get_orthonormal_mapping(mapvars, cam, ah);
    
    mapvars = map_wide2pers_normal_bestfit(img, mapvars, ah);
    mapvars = fisheye_clickborder(cam, mapvars, ah, false);
    mapvars = get_orthonormal_mapping(mapvars, cam, ah);
    mapvars = map_wide2pers_normal_bestfit(img, mapvars, ah);
    gpoints = fisheye_clickborder(cam, mapvars, ah, true);
   
    % Find grid points
    gpoints = auto_grid_find(mapvars, gpoints, ah);
        
    % The grid points are for the normalised and rescaled image.  Need to
    % remap back the the original image plane coords (u,v).
    gpoints_original = remap_grid(gpoints, mapvars);
    
    % May need to transpose the grid order 
    % (just the convention used for grid detection)
    %[nr,nc] = size(gpoints_original(:,:,1))
    %if nr ~= (numgrid(2)+1)
    gpoints_original = permute(gpoints_original, [2 1 3]);
    %end
    
    % Rough initial estimate of extrinsic rotation
    Rnormals = mapvars.normal;
end


%
% Setup figure
% Spanning subplots to keep things simple
%
function ah = init_grid_detect_figure(nr, nc)
  
    fh = figure(1);
    colormap(gray(255))
    %ah(1) = subplot(4,6,[1:4 7:10 13:16 19:22]);
    ah(1) = subplot(5,7,[1:5 8:12 15:19 22:26 29:33]);
    axis([0 nc  0 nr]);
    %axis off
    ah(2) = subplot(5,7,[6 7]);
    ah(3) = subplot(5,7,[13 14 20 21]);
    ah(4) = subplot(5,7,[27 28 34 35]);
    %ah(2) = subplot(4,6,[5 6]);
    %ah(3) = subplot(4,6,[17 18]);
    %ah(4) = subplot(4,6,[23 24]);
    %axis off
    set(fh, 'Position', [100 300 1300 1300/2.2], 'Children', ah)
    title(ah(2),'Greyscale Histogram','FontWeight','Bold')     
    drawnow
end

