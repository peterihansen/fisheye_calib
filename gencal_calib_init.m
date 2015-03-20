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
    %================================================================
    state = zeros(length(igrid)*7, 1);
    for i = 1:length(igrid)
        ind = (i-1)*7 + 1;
        state(ind:(ind+6)) = [calib.q(:,i); calib.t(:,i)];
    end

    
    % The optimization
    %================================================================
    options = optimset('lsqnonlin');
    %options.OutPutFcn = '
    options.Algorithm = 'levenberg-marquardt';
    options.Display = 'iter';
    options.MaxFunEvals = 5e4;
    options.MaxIter = 5e4;
    options.JacobPattern = jacobpattern_extrinsic(igrid);
    state = lsqnonlin(@ofun_gencal_calib_init,state,[],[],options);
    
    
    
    % The objective function
    %================================================================
    function err = ofun_gencal_calib_init(state)
       
        % Preallocate the matrices A and b
        %--------------------------------------------------
        n = tgrid.nx * tgrid.ny;
        A = zeros(2*n*length(igrid), calib.porder+2);
        b = zeros(2*n*length(igrid), 1);
        
        
        % For each image, get the spherical polar angles
        % and add to the matrix A
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
                A(ind_rows_x,k) = (cos(phi).*(theta.^k))';
                A(ind_rows_y,k) = (sin(phi).*(theta.^k))';
            end
            A(ind_rows_x,k+1) = ones(n,1);
            A(ind_rows_y,k+2) = ones(n,1);
            
            % Update b
            b(ind_rows_x,1) = igrid(j).u';
            b(ind_rows_y,1) = igrid(j).v';
        end
        
        % Solve the linear system and get the errors
        %--------------------------------------------------
        %A = [A(:,1) A(:,3:end)];
        poly_pp = A\b;
        err = A*poly_pp - b;
        %poly_pp = [poly_pp(1); 0; poly_pp(2:end)];
        
%         figure(2);cla;plot(err);drawnow
%         fprintf('Error:  mean abs = %0.4f  std = %0.4f\n', ...
%                     sqrt(sum(err.^2)) / (n*length(igrid)), std(err));
    end
    %================================================================
    
    
    
    
    % Store result for extrinsic and update calib
    %================================================================
    for i = 1:length(igrid)
        ind = (i-1)*7 + 1;
        q = state(ind:(ind+3));
        t = state((ind+4):(ind+6));
        calib.q(:,i) = q / sqrt(sum(q.^2));
        calib.t(:,i) = t;
    end
    calib.pvars = poly_pp(1:calib.porder);
    pp = poly_pp((calib.porder+1):end);
    calib.H = [1 0 pp(1); 0 1 pp(2); 0 0 1];
end



function jacobpattern = jacobpattern_extrinsic(igrid)

    N = length(igrid);
    n = length(igrid(1).u);

    rows = zeros(7*2*N*n,1);
    cols = zeros(7*2*N*n,1);
    
    [col_base, row_base] = meshgrid(1:7,1:(2*n));
    col_base = col_base(:);
    row_base = row_base(:);
    ind_base = 1:(7*2*n);
    
    for i = 1:N
        ind_row = (i-1)*2*n;
        ind_col = (i-1)*7;
        rows(ind_base,1) = row_base + ind_row;
        cols(ind_base,1) = col_base + ind_col;
        ind_base = ind_base + 7*2*n;
    end
    jacobpattern = sparse(rows, cols, ones(7*2*N*n,1));
end





