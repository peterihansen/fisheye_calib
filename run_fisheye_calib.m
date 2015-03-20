%
%   Run generic polynomial calibration
%
%   Works with wide-angle fisheye
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


clear all
% addpath(genpath('./'));
% addpath(genpath('../'));

% OPTIONS
%====================================================================

% Root directory
%datadir = '/media/foxdata1/data/vslam/2013_05_20_EC/calibration/fisheye_intrinsic';

datadir = '~/proj/vslam/data/2014_09_01_EC/calibration/fisheye_intrinsic';
imgext = 'pgm';
isbayer = true;

% Number of grids (x,y) and spacing ratio
numgrid = [9; 7];
ratio_YonX = 1.0;

% Setup an initial approx camera calibration and 
% include function handles for the mappings
cam_init.u = 600;
cam_init.v = 595;
cam_init.r90 = 500;
cam_init.map_img2ray = @map_unified_img2sphere;
cam_init.map_X2img = @map_unified_X2img;

% Polynomial order
porder = 5;

%====================================================================



% Get approx unified image model params which will be used
% to produce rectified orthonormal views during grid detection
%----------------------------------------------------------
cam_init = set_unified_model(cam_init);


% STEP 1: Grid detection
% Keep prompting for new images, load, detect grid points, and save
%====================================================================
imgdir = strcat(datadir,'/images');
griddir = strcat(datadir,'/gridpoints');
loop_grid = true;
while loop_grid == true;
    
    % Prompt for new image
    %-----------------------------------------------------------
    reply = input('Include another image ? y/n [y]: ', 's');
    if ~isempty(reply) & reply ~= 'y'
        loop_grid = false;
        continue
    end
    
    % Load the image and convert greyscale (assuming 8bit)
    %-----------------------------------------------------------
    img_name = uigetfile(sprintf('%s',imgdir,'/*.*'), 'Select Image');
    img_file = sprintf('%s',imgdir,'/',img_name);
    [pathtmp,img_name,ext] = fileparts(img_name);
    img = imread(img_file);
    if (ndims(img) == 3)
        img = rgb2gray(img);
    elseif isbayer == true
        img = rgb2gray(demosaic(img,'rggb'));
    end
    [nr nc] = size(img);
    img = double(img);
    
    % User guided grid detection
    % Returns grid points and estimated extrinsic pose
    %-----------------------------------------------------------
    [gpoints, Rext] = image_grid_detect(img, cam_init, numgrid);
    
    figure(2)
    colormap(gray(255))
    cla; imagesc(img); axis image; hold on
    plot(gpoints(:,:,1)+1, gpoints(:,:,2)+1, 'yo','Markerfacecolor','r','Markersize',5);

    % Save the grid points
    %-----------------------------------------------------------
    svfile = strcat(griddir,sprintf('/gpoints.%s.mat',img_name));
    save(svfile,'gpoints', 'Rext');
end



% STEP 2: Reload all grid points, and setup the scene points
% Store them in vector formats
%====================================================================
gpoints_fnames = glob(strcat(griddir,'/gpoints*.mat'));
Nimg = length(gpoints_fnames);
Rext_init = zeros(3,Nimg);
for i = 1:Nimg
    load(gpoints_fnames{i});
    utmp = gpoints(:,:,1);
    vtmp = gpoints(:,:,2);
    igrid(i).u = utmp(:)';
    igrid(i).v = vtmp(:)';
    Rext_init(:,i) = Rext;
end
tgrid.nx = numgrid(1)+1;
tgrid.ny = numgrid(2)+1;
tgrid.dx = 0.1;
tgrid.dy = ratio_YonX * tgrid.dx;
[X,Y] = meshgrid((0:(tgrid.nx-1))*tgrid.dx, (0:(tgrid.ny-1))*tgrid.dy);
tgrid.X = [X(:)';Y(:)';zeros(1,tgrid.nx*tgrid.ny)];


% STEP 3: Do the calibration
%====================================================================
cam = initial_checkerboard_poses(cam_init, tgrid, igrid, Rext_init);

cam.porder = porder;
cam = gencal_calib_init(cam, tgrid, igrid);
cam = gencal_calib_refine(cam, tgrid, igrid);
%cam = rmfield(cam,'poly');

% Add function handles for the mappings
cam.map_img2ray = @polyH_map_img2sphere;
cam.map_ray2img = @polyH_map_X2img;


% Look at the results/errors
%====================================================================
analyse_calibration(cam, tgrid, igrid);


% Save the result as matlab struct
%====================================================================
save(strcat(datadir,'/cam.mat'), 'cam','tgrid','igrid');




