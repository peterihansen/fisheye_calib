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

function morpho = find_intersection(patch)

    [nr,nc] = size(patch);
    morpho.patch = patch;
    
    % Histogram the values to get the binary threshold value
    %[h] = ihist(patch); 
    h = 1;
    thresh = mean(mean(patch));
    morpho.h = h;
    morpho.thresh = thresh;
    
    % Threshold the greyscae values
    bpatch = ones(nr,nc);
    [i,j] = find(patch < thresh);
    for k = 1:length(i)
        bpatch(i(k),j(k)) = 0;
    end
    
    % Apply morphological operator to remove some edge noise
    for k = 1:3
        temp_patch = bpatch;
        bpatch(2:nr-1,2:nc-1) = 0;
        for i=2:nr-1
            for j=2:nc-1
                z = find(temp_patch(i-1:i+1,j-1:j+1) == 1);
                if length(z) >= 4
                    bpatch(i,j) = 1;
                else
                    bpatch(i,j) = 0;
                end
            end
        end
    end
    
    %figure(66);cla;imagesc(bpatch);pause
    
    % Label the blobs in the binary patch
    lpatch = ilabel(bpatch);
    morpho.lpatch = lpatch;    
    nblobs = max(max(lpatch));
    
    
    % Get the position of the feature
    %====================================================================
    
    % LESS THAN 3 REGIONS OR MORE THAN 4 REGIONS (return error)
    if ((nblobs <= 2) | (nblobs > 4))
        morpho.flag = 1;
        return
    
    % EXACTLY FOUR REGIONS - FEATURE AT THE INTERSECTION OF ALL    
    elseif nblobs == 4
        morpho.flag = 0;
        for i = 1:nr-1
            for j = 1:nc-1
                h = hist([lpatch(i,j:j+1) lpatch(i+1,j:j+1)], [1 2 3 4]);
                if h == [1 1 1 1]; 
                    morpho.u = i + 0.5
                    morpho.v = j + 0.5
                    morpho.p1 = [morpho.u morpho.v];
                    morpho.p2 = morpho.p1;
                    return
                end
            end
        end
        morpho.flag = 1;
    
        
    % THREE REGIONS - FIND MIN DISTANCE BETWEEN TWO NON-SEPARATING REGIONS
    else 
        morpho.flag = 0; 
        
        % Put the border values into a single vector
        lborder = [lpatch(2,1:nc) lpatch(2:nr,nc)' ...
                   lpatch(nr:nc-1:-1:1) lpatch(nr-1:-1:1,1)'];
   
        % Ensure there is no wrap around value
        border_length = length(lborder); 
        if lborder(1) == lborder(border_length)
            val = lborder(1);
            index = find(lborder ~= val, 1, 'first');
            lborder_temp = lborder;
            lborder = [lborder_temp(index:border_length), ...
                       lborder_temp(1:index-1)];
        end
    
        % 'Blob' the border vector and find the separating region
        border_blob(1) = lborder(1);
        index_old = 1;
        for i = 1:3
            index = find(lborder(index_old:border_length) ~= ...
                                                border_blob(i),1,'first');
            index_old = index + index_old - 1;
            border_blob(i+1) = lborder(index_old);
        end
        h = hist(border_blob, [1 2 3]);
        [val,sep_val] = max(h);
    
    
        % Get the valid border points of the non-separating regions
        vals = find([1 2 3] ~= sep_val);
        k1 = 0;
        k2 = 0;
        for i = 2:nr-1
            vect = lpatch(i,2:nc-1);
            x1f = find(vect==vals(1), 1, 'first');
            x1l = find(vect==vals(1), 1, 'last');
            x2f = find(vect==vals(2), 1, 'first');
            x2l = find(vect==vals(2), 1, 'last');
        
            % Region 1
            if isempty(x1f)~=1
                for j = x1f+1:x1l+1
                    if ((lpatch(i,j-1)~=vals(1)) | ...
                        (lpatch(i,j+1)~=vals(1)) | ...
                        (lpatch(i-1,j)~=vals(1)) | ...
                        (lpatch(i+1,j)~=vals(1)))
                        k1 = k1 + 1;
                        blob_border1(k1,:) = [j i];
                    end
                end
            end

            % Region 2
            if isempty(x2f)~=1
                for j = x2f+1:x2l+1
                    if ((lpatch(i,j-1)~=vals(2)) | ...
                        (lpatch(i,j+1)~=vals(2)) | ...
                        (lpatch(i-1,j)~=vals(2)) | ...
                        (lpatch(i+1,j)~=vals(2)))
                        k2 = k2 + 1;
                        blob_border2(k2,:) = [j i];
                    end
                end
            end       
        end
        morpho.blob_border1 = blob_border1;
        morpho.blob_border2 = blob_border2;
    
    
        % Find the minimun separation between the blobs
        for i = 1:length(blob_border1(:,1))
            for j = 1:length(blob_border2(:,2))
                dist(i,j) = sqrt(sum((blob_border1(i,:) - ...
                                            blob_border2(j,:)).^2));
            end
        end
    
        % Find the border points on each region with the smallest distance
        min_dist = min(min(dist));
        [index_i, index_j] = find(dist == min_dist);
    
        % Check that there is not multiple minima
        if length(index_i) > 1
            for i = 1:length(index_i)
                p1(i,:) = blob_border1(index_i(i),:);
                p2(i,:) = blob_border2(index_j(i),:);
                u(i) = (p1(i,1) + p2(i,1)) / 2;
                v(i) = (p1(i,2) + p2(i,2)) / 2;
            end
            morpho.u = mean(u);
            morpho.v = mean(v);
        else
            p1 = blob_border1(index_i,:);
            p2 = blob_border2(index_j,:);
            morpho.u = (p1(1,1) + p2(1,1)) / 2;
            morpho.v = (p1(1,2) + p2(1,2)) / 2;
        end
        morpho.p1 = p1;
        morpho.p2 = p2;
    end
    