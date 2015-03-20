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

function gpoints = auto_grid_find(imvars, gpoints_temp, ah)

    axes(ah(1))
    title(ah(1),'Finding Intersections','FontWeight','Bold', 'FontSize',12)
    
    % Get image size
    [nr,nc] = size(imvars.remap_img);

    % Make sure gpoints is in the correct order (min-max)
    [gr,gc] = size(gpoints_temp(:,:,1));
    [val,index] = max([gr gc]);    
    if index == 2
        gpoints(:,:,1) = gpoints_temp(:,:,1)';
        gpoints(:,:,2) = gpoints_temp(:,:,2)';
    else
        gpoints = gpoints_temp;
    end
    
    
%     % Save the points and remapped image for plotting later
%     svfile = 'points_est.mat';
%     save(svfile,'gpoints');
%     svfile = 'img_remap2.pgm';
%     imwrite(imvars.remap_img/255,svfile,'pgm');
    
    
    % Try finding all the intersections automatically
    missing = 0;
    imax = length(gpoints(:,1,1));
    jmax = length(gpoints(1,:,1));
    for i = 1:imax;
        for j = 1:jmax
        
            U = round(gpoints(i,j,1));
            V = round(gpoints(i,j,2));
            Pw = find_patch_widths(gpoints, i, imax, j, jmax, nr, nc); 
            patch = imvars.remap_img(V-Pw(3):V+Pw(4), U-Pw(1):U+Pw(2));
            
            
            morpho = find_intersection(patch);
        
            % Plot the binary thresholding histogram
            axes(ah(2))
            cla
            hold on
            h = morpho.h / sum(morpho.h);
            bar(h);
            plot([morpho.thresh morpho.thresh],[0 max(h)],'r')
            axis([-1 256 0 max(h)])
        
            % Plot the input patch
            axes(ah(3))
            cla
            image(patch)
            axis image
            hold on
            plot(morpho.u, morpho.v, 'y.', 'markersize',15);
        
            % Plot the binary patch
            axes(ah(4))
            cla
            image(255 * morpho.lpatch / max(max(morpho.lpatch)))
            axis image
            hold on
            if isfield(morpho,'blob_border1')
                for k = 1:length(morpho.blob_border1)
                    plot(morpho.blob_border1(k,1), ...
                        morpho.blob_border1(k,2),'g.','linewidth',2)
                end
            end
            if isfield(morpho,'blob_border2')
                for k = 1:length(morpho.blob_border2)
                    plot(morpho.blob_border2(k,1), ...
                        morpho.blob_border2(k,2),'b.','linewidth',2)
                end
            end
    
            for k = 1:length(morpho.p1(:,1))
                plot([morpho.p1(k,1) morpho.p2(k,1)], ...
                     [morpho.p1(k,2) morpho.p2(k,2)], 'r', 'linewidth',2)
            end
            
        
            % Update the points 
            if morpho.flag == 0
                gpoints(i,j,1) = morpho.u + min([U-Pw(1),U+Pw(2)]) - 1;
                gpoints(i,j,2) = morpho.v + min([V-Pw(3),V+Pw(4)]) - 1;
                gpoints = update_grid(gpoints, i, imax, j, jmax);
            else
                missing = 1;
                gpoints(i,j,1) = NaN;
                gpoints(i,j,2) = NaN;
            end
                
            axes(ah(1))
            plot(gpoints(i,j,1),gpoints(i,j,2),'y.','markersize',15)
            drawnow
        end
    end
    
    
    % 'Fill in' any NaN values (features that were not found)
    if missing == 1
        [nr, nc] = size(gpoints(:,:,1));
        [i,j] = find(isnan(gpoints(:,:,1)));
        
        for k = 1:length(i)    
            a = 0;
            b = 0;
            clear xr xc yr yc
            
            % Get the row points
            c = j(k);
            for r = i(k)-1:2:i(k)+1
                if ((r >= 1) & (r <= nr))
                    a = a + 1;
                    xr(a) = gpoints(r,c,1);
                    yr(a) = gpoints(r,c,2);
                end
            end
            
            % Get the column points
            r = i(k);
            for c = j(k)-1:2:j(k)+1
                if ((c >= 1) & (c <= nc))
                    b = b + 1;
                    xc(b) = gpoints(r,c,1);
                    yc(b) = gpoints(r,c,2);
                end
            end
            
            % Get the estimates of the midpoints
            if length(xr) == 2
                midxr = sum(xr) / 2;
                midyr = sum(yr) / 2;
                addr = 1;
            else
                midxr = 0;
                midyr = 0;
                addr = 0;
            end
            
            if length(xc) == 2
                midxc = sum(xc) / 2;
                midyc = sum(yc) / 2;
                addc = 1;
            else
                midxc = 0;
                midyc = 0;
                addc = 0;
            end
            
            u = (midxr*addr + midxc*addc) / (addr + addc);
            v = (midyr*addr + midyc*addc) / (addr + addc);
            
            gpoints(i(k),j(k),1) = u;
            gpoints(i(k),j(k),2) = v;
            
            plot(u,v,'bp')
        end
    end
end
    