% Get the user to click on the grid corners
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

function mapvars = checkerboard_clickborder(cam, mapvars, ah, flagdone)
    
    axes(ah(1))
    cla
    image(mapvars.remap_img)
    axis image
    hold on
    
    % Click the 4 outermost corners
    title(ah(1),'Click 4 Outermost Corners: clockwise, x direction first',...
                                    'FontWeight','Bold', 'FontSize',12)
    
    % The first time, simply pick near the grid border corners
    if mapvars.remapped == 0
        for i = 1:4
            [u(i),v(i)] = ginput(1);
            plot(u(i),v(i),'rx','linewidth',2)
        end
        points = [u' v']-1;
        mapvars.border_wide = parallel_lines(points);
        
    else
        
        % Reproject the sets of lines found onto the remapped image
        for k = 1:2
            for j = 1:2
                
                S = cam.map_img2ray(cam, mapvars.border_wide{k}{j}');
                Sn = map_sphere2sphere(S, mapvars.normal, 1);
                U = mapvars.remap_K * Sn;
                U = U(1:2,:) ./ repmat(U(3,:),[2 1]);
                plines{k}{j} = U';
                
%                 for i = 1:length(mapvars.border_wide{k}{j})
%                     
%                     
%                     x = mapvars.border_wide{k}{j}(i,1) - cam.Cx;
%                     y = mapvars.border_wide{k}{j}(i,2) - cam.Cy;
%                     S = map_unified_img2sphere(x, y, cam.m, cam.l); 
%                     Sn = map_sphere2sphere(S, mapvars.normal, 1);
%                     
%                     % PINHOLE BIT
%                     
%                     [X,Y] = map_unified_sphere2img(Sn(1), Sn(2), ...
%                                                 Sn(3), cam.m, 0);
%                                             
%                     % Correct for the origin and scaling
%                     U = (X - mapvars.remap_originX) / mapvars.remap_spacing;
%                     V = (Y - mapvars.remap_originY) / mapvars.remap_spacing;
%                     
%                     plines{k}{j}(i,1) = U;
%                     plines{k}{j}(i,2) = V;
%                 end
            end
        end
            
        % Get the estimated grid coords and overlay
        gpoints = grid_coords(plines, mapvars.numgrid, 0);
        
        if flagdone == true
            for i = 1:length(gpoints(:,1,1))
                for j = 1:length(gpoints(1,:,1))
                    plot(gpoints(i,j,1)+1,gpoints(i,j,2)+1,'rx', 'linewidth',2)
                end
            end
            %clear reply
            %reply = input('Use this estimate y/n [y]?: ','s');
            %if reply == 'y'
            mapvars = gpoints;
            return
        end
        
        
        % Select the new points using the assisted feature selection
        cla
        image(mapvars.remap_img)
        axis image
        hold on
        title(ah(1),'Reselect the points: : clockwise, x direction first',...
                                    'FontWeight','Bold', 'FontSize',12)
               
        i = 0;
        U = zeros(2,4);
        
        while i < 4  %8
            axes(ah(1))
            [Uc,Vc] = ginput(1);
            
            % Take a small patch and get position of feature in patch
            %patch = mapvars.remap_img(round(Vc(1)):round(Vc(2)), ...
            %                          round(Uc(1)):round(Uc(2)));
                               
            % Find the intersection automatically
            %morpho = find_intersection(patch);
            
            %if morpho.flag == 0
            i = i + 1;
                
                % Convert back to the image coordinates
            %    U(i) = morpho.u + round(Uc(1)) - 1 - 1;
            %    V(i) = morpho.v + round(Vc(1)) - 1 - 1;
            U(:,i) = [Uc;Vc];
            
            % Plot the position on the image
            plot(U(1,i), U(2,i),'rx','linewidth',2)
        
            % Plot the automatic intersection info (morphology)
            %plot_morpho(morpho, ah);
                       
%             % Correct for any offset of origin or rescaling
%             X = U(i) * mapvars.remap_spacing + mapvars.remap_originX;
%             Y = V(i) * mapvars.remap_spacing + mapvars.remap_originY;
%                 
%             % Map point to sphere
%             Sn = map_unified_img2sphere(X, Y, mapvars.m, 0);
%               
%             % De-correct normal
%             S(i,:) = map_sphere2sphere(Sn, mapvars.normal, 0);
%                 
%             % Map to wide angle image
%             [x(i) y(i)] = map_unified_sphere2img(S(i,1), S(i,2), ...
%                                             S(i,3), mapvars.m, mapvars.l);
%             u(i) = x(i) + mapvars.Cx;
%             v(i) = y(i) + mapvars.Cy;
%             %end
        end

        % Map back to original image
        S = inv(mapvars.remap_K) * [U; ones(1,4)];
        S = S ./ repmat(sqrt(sum(S.^2)),[3 1]);
        S = map_sphere2sphere(S, mapvars.normal, 0);
        U = cam.map_ray2img(cam, S);
        points = U';       %[u' v'];
        
        % Get the sets of parallel lines
        mapvars.border_wide = parallel_lines(points);
        output = mapvars;
    end
    drawnow 
end


%         % Select the new points using the assisted feature selection
%         cla
%         image(imvars.remap_img)
%         axis image
%         hold on
%         title(ah(1),'Click 8 Points on Border in Order',...
%                                     'FontWeight','Bold', 'FontSize',12)
%                
%         i = 0;
%         while i < 8
%             axes(ah(1))
%             [Uc,Vc] = ginput(2);
%             
%             % Take a small patch and get position of feature in patch
%             patch = imvars.remap_img(round(Vc(1)):round(Vc(2)), ...
%                                      round(Uc(1)):round(Uc(2)));
%                                
%             % Find the intersection automatically
%             morpho = find_intersection(patch);
%             
%             if morpho.flag == 0
%                 i = i + 1;
%                 
%                 % Convert back to the image coordinates
%                 U(i) = morpho.u + round(Uc(1)) - 1 - 1;
%                 V(i) = morpho.v + round(Vc(1)) - 1 - 1;
%                 
%                 % Plot the position on the image
%                 plot(U(i)+1,V(i)+1,'ro','linewidth',2)
%         
%                 % Plot the automatic intersection info (morphology)
%                 plot_morpho(morpho, ah);
%                        
%                 % Correct for any offset of origin or rescaling
%                 X = U(i) * imvars.remap_spacing + imvars.remap_originX;
%                 Y = V(i) * imvars.remap_spacing + imvars.remap_originY;
%                 
%                 % Map point to sphere
%                 Sn = map_unified_img2sphere(X, Y, imvars.m, 0);
%               
%                 % De-correct normal
%                 S(i,:) = map_sphere2sphere(Sn, imvars.normal, 0);
%                 
%                 % Map to wide angle image
%                 [x(i) y(i)] = map_unified_sphere2img(S(i,1), S(i,2), ...
%                                             S(i,3), imvars.m, imvars.l);
%                 u(i) = x(i) + imvars.Cx;
%                 v(i) = y(i) + imvars.Cy;
%             end
%         end
% 
%         points = [u' v'];
%         
%         % Get the sets of parallel lines
%         imvars.border_wide = parallel_lines(points) ;
%         output = imvars;
%     end
%     
%     drawnow
    