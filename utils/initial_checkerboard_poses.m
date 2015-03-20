%
%   Get the pose P for each of the checkerboard images
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

function calib = initial_checkerboard_poses(cam, tgrid, igrid, Rext)

    Xg = tgrid.X;
    npoints = length(igrid(1).u);
    
    for i = 1:length(igrid)
    
        % Spherical coordinates of the grid points
        U = [igrid(i).u; igrid(i).v];
        S = cam.map_img2ray(cam, U);
        
        %for j = 1:length(igrid(i).u)
            %S(:,j) = map_poly_homog_img2sphere(U, cam.H, cam.pvars);
        %    S(:,j) = map_unified_img2sphere(U(1,j)-cam.Cx,U(2,j)-cam.Cy, cam.m, cam.l);
        %end
        
        ea = Rext(:,i);
        %R = euler_matrix('z',ea(3)) * euler_matrix('y',ea(2)) * ...
        %    euler_matrix('z',ea(1) + pi/2);
        %-------------------------------------------------
        beta = Rext(1,i);
        gamma = Rext(2,i);
        alpha = Rext(3,i);
        
        loop = true;
        while loop == true
            Ry = euler_matrix('y',beta);
            Rz = euler_matrix('z',gamma);
            Rz2 = euler_matrix('z',alpha);
            R = Rz * Ry * Rz2;
            %--------------------------------------------------
        
            % Rotate spherical and proj to plane
            Stmp = R' * S;
            Stmp = Stmp ./ repmat(Stmp(3,:),[3 1]);
            
            % The translation
            tz =  (sqrt(sum((Xg(1:2,end) - Xg(1:2,1)).^2)) / ...
                   sqrt(sum((Stmp(1:2,end) - Stmp(1:2,1)).^2)));
            Stmp = tz * Stmp;
            tx = mean(Xg(1,:) - Stmp(1,:));
            ty = mean(Xg(2,:) - Stmp(2,:));
            Stmp(1,:) = Stmp(1,:) + tx;
            Stmp(2,:) = Stmp(2,:) + ty;
            dl = sqrt(sum((Stmp(1:2,1) - Xg(1:2,1)).^2));
            if dl > (2*abs(Xg(2,1) - Xg(2,2)))
                alpha = alpha + pi;
            else
                loop = false;       % Remove this loop - could be infinite
            end
        end
        
%         figure(10)
%         cla;hold on
%         plot(Xg(1,:),Xg(2,:),'b.')
%         plot(Stmp(1,:),Stmp(2,:),'r.')
%         for n = 1:length(Stmp(1,:))
%             plot([Xg(1,n) Stmp(1,n)], [Xg(2,n) Stmp(2,n)],'-g')
%         end
%         axis equal

        % The final extrinsic camera pose
        P = [R -R*([tx;ty;-tz]); 0 0 0 1];
        q = rot2quat(P(1:3,1:3));
        t = P(1:3,4);
        calib.q(:,i) = q;
        calib.t(:,i) = t;
            
%         % Test
%         Xcam = P * [Xg;ones(1,length(Xg(1,:)))];
%         Scam = Xcam(1:3,:) ./ repmat(sqrt(sum(Xcam(1:3,:).^2)),[3 1]);
%         U2 = map_unified_sphere2img(Scam(1,:), Scam(2,:), Scam(3,:), cam.m, cam.l);
%         U2(1,:) = U2(1,:) + cam.Cx;
%         U2(2,:) = U2(2,:) + cam.Cy;
        
        
%         figure(11);
%         cla; hold on
%         plot3(S(1,:),S(2,:),S(3,:),'b.')
%         plot3(Scam(1,:),Scam(2,:),Scam(3,:),'r.')
%         
%         
%         
%         figure(12);
%         cla;hold on
%         plot(U(1,:),U(2,:),'b.')
%         plot(U2(1,:),U2(2,:),'r.')
%         for n = 1:length(U(1,:))
%            plot([U(1,n) U2(1,n)],[U(2,n) U2(2,n)],'-g')
%         end
%             
%         
%         pause
    end
end



function U = map_unified_sphere2img(x, y, z, m, l)
    scale_factor = (l+m) ./ (l+z);
    U(1,:) = x .* scale_factor; 
    U(2,:) = y .* scale_factor;
end



