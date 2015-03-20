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

function [calib,igrid] = gencal_calib_init(calib, tgrid, igrid)


    % The state vector
    % [extrinsic; matrix H; poly]
    %================================================================
    
    % Extrinsic
    state = zeros(length(igrid)*7, 1);
    for i = 1:length(igrid)
        ind = (i-1)*7 + 1;
        state(ind:(ind+6)) = [calib.q(:,i); calib.t(:,i)];
    end
    
    % Camera H
    state = [state; 1; 0];
    
    
    % The optimization
    %================================================================
    options = optimset('lsqnonlin');
    options.Algorithm = 'levenberg-marquardt';
    options.Display = 'iter';
    options.MaxFunEvals = 1e4;
    options.MaxIter = 1e4;
    state = lsqnonlin(@ofun_gencal_calib_refine,state,[],[],options);
    
    
    
    % The objective function
    %================================================================
    function err = ofun_gencal_calib_refine(state)
       
        % Preallocate the matrices A and b
        %--------------------------------------------------
        n = tgrid.nx * tgrid.ny;
        A = zeros(2*n*length(igrid), calib.porder + 2);
        b = zeros(2*n*length(igrid), 1);
       
        
        % Camera matrix affine component
        ind = length(igrid)*7 + 1;
        H_y = state(ind);
        H_s = state(ind+1);
             
        % For each image, project points to image and
        % get the error
        %--------------------------------------------------
        for j = 1:length(igrid)
            ind = (j-1)*7 + 1;
            q = state(ind:(ind+3));
            R = quat2rot(q / sqrt(sum(q.^2)));
            t = state((ind+4):(ind+6));
            Xp = R * tgrid.X + repmat(t, [1 n]);
            Xp = Xp ./ repmat(sqrt(sum(Xp.^2)), [3 1]);
            theta = acos(Xp(3,:));
            phi = atan2(Xp(2,:),Xp(1,:));
        
            ind_A = (j-1)*2*n + 1;
            ind_rows_x = ind_A:2:(ind_A+2*n-1);
            ind_rows_y = (ind_A+1):2:(ind_A+2*n-1);
           
            % Update A
            for k = 1:calib.porder
                %A(ind_rows_x,k) = (cos(phi).*(theta.^k))';
                %A(ind_rows_y,k) = (sin(phi).*(theta.^k))';
                A(ind_rows_x,k) = ((theta.^k) .* (cos(phi) + sin(phi)*H_s))';
                A(ind_rows_y,k) = ((theta.^k) .* (sin(phi)*H_y))';
            end
            A(ind_rows_x,k+1) = ones(n,1);
            A(ind_rows_y,k+2) = ones(n,1);
            
            % Update b
            b(ind_rows_x,1) = igrid(j).u';
            b(ind_rows_y,1) = igrid(j).v';
        end
        
        
        % Solve the linear system and get the errors
        %--------------------------------------------------
        poly_pp = A\b;
        err = A*poly_pp - b;
        
        %fprintf('Error:  mean abs = %0.8f  std = %0.4f\n', ...
        %                sqrt(sum(err.^2)) / (2*length(err)), std(err));
       
    end
    %================================================================
    
    
    
    % Store result for extrinsic and update calib
    %================================================================
    
    % Extrinsic
    for i = 1:length(igrid)
        ind = (i-1)*7 + 1;
        q = state(ind:(ind+3));
        t = state((ind+4):(ind+6));
        calib.q(:,i) = q / sqrt(sum(q.^2));
        calib.t(:,i) = t;
    end
    
    % Polynomial
    calib.poly = poly_pp(1:calib.porder);
    
    % Camera matrix H
    ind = length(igrid)*7 + 1;
    Hvars = state(ind:(ind + 1));
    calib.H = [1 Hvars(2) poly_pp(end-1); 0 Hvars(1) poly_pp(end); 0 0 1];
    %calib.H = [Hvars' poly_pp(end-1); Hvars(2) 1 poly_pp(end); 0 0 1];      
end





