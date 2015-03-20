% 
%   Show calibration results - errors etc.
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

function analyse_calibration(cam, tgrid, igrid)

    N = length(igrid);
    n = tgrid.nx * tgrid.ny;

    % Get the errors for each image
    %======================================================
    err_u = zeros(N, n);
    err_v = zeros(N, n);
    rad = zeros(N,n);
    for i = 1:N
        q = cam.q(:,i);
        R = quat2rot(q);
        t = cam.t(:,i);
        Sp = R * tgrid.X + repmat(t, [1 n]);
        Sp = Sp ./ repmat(sqrt(sum(Sp.^2)), [3 1]);
        Up = cam.map_ray2img(cam, Sp);
        err_u(i,:) = igrid(i).u - Up(1,:);
        err_v(i,:) = igrid(i).v - Up(2,:);
        rad(i,:) = sqrt((igrid(i).u - cam.H(1,3)).^2 + ...
                        (igrid(i).v - cam.H(2,3)).^2);
    end
    sum(err_u(:).^2 + err_v(:).^2)
    
    
    % Plot all together: colorize different images
    %===========================================================
    figure(1)
    clf; hold on
    cmap = colormap(lines(N));
    for i = 1:N
        plot(err_u(i,:), err_v(i,:), 'bx', 'Color', cmap(i,:));
    end
    plot(mean(err_u(:)), mean(err_v(:)), 'k.','MarkerSize',20)
    alim = ceil(2*max(abs([err_u(:); err_v(:)])))/2;
    set(gca,'XLim',[-alim alim],'YLim',[-alim alim]);
    set(gca,'XTick',-alim:0.5:alim,'YTick',-alim:0.5:alim);
    set(gca,'FontSize',14)
    title(sprintf('Mean = (%0.2f , %0.2f)', mean(err_u(:)), mean(err_v(:))),'FontSize',14)
    xlabel('x error (pixels)','FontSize',14)
    ylabel('y error (pixels)','FontSize',14)
    set(gcf,'PaperPositionMode','Auto')
    grid on; box on
    
    
    % For each image, plot mean u,v errors and std. u,v errors
    % Use boxplot (like the stereo extrinsic calibration)
    %===========================================================
    figure(3)
    cla; hold on
    subplot(2,1,1)
    errorbar(mean(err_u,2), 2*std(err_u,[],2), 'bx')
    xlabel('Image Number','FontSize',16)
    ylabel('Value (pixels)','FontSize',16)
    title('Mean u error with 2\sigma error bars','FontSize',16)
    set(gca,'FontSize',16)
    grid on
    subplot(2,1,2)
    errorbar(mean(err_v,2), 2*std(err_v,[],2), 'bx')
    xlabel('Image Number','FontSize',16)
    ylabel('Value (pixels)','FontSize',16)
    title('Mean v error with 2\sigma error bars','FontSize',16)
    set(gca,'FontSize',16)
    grid on
    drawnow
    
%     figure(2)
%     clf;
%     subplot(1,2,1);cla;hold on
%     plot(mean(err_u,2), '-b.')
%     plot(mean(err_v,2), '-r.')
%     legend('u','v')
%     set(gca,'FontSize',14)
%     title(sprintf('Mean u,v pixels errors'))
%     xlabel('Image Number','FontSize',14)
%     ylabel('Error (pixels)','FontSize',14)
%     grid on
%     subplot(1,2,2);cla;hold on
%     plot(std(err_u,[],2),'-b.')
%     plot(std(err_v,[],2),'-r.')
%     legend('u','v')
%     set(gca,'FontSize',14)
%     title(sprintf('Standard deviation of u,v pixels errors'))
%     xlabel('Image Number','FontSize',14)
%     ylabel('Standard Deviation (pixels)','FontSize',14)
%     grid on
%     set(gcf,'PaperPositionMode','Auto')
    
    




    % Plot pdf of the errors (u,v,all)
    %===========================================================
    ex = err_u(:);   % Can put all together now
    ey = err_v(:);
    
    alim = max(abs([ex;ey]));
    alim = ceil(alim * 5) / 5
    xbins = -alim:0.25:alim;
    ybins = -alim:0.25:alim;
    hx = hist(ex, xbins);
    hy = hist(ey, ybins);
    
    figure(10);
    cla;
    hold on
    bar(xbins, hx/sum(hx));
    x_mean = mean(ex);
    x_std = std(ex);
    xg = -alim:0.1:alim;
    g = exp(-(xg-x_mean).^2 / (2*x_std^2));
    gtmp = exp(-(xbins-x_mean).^2 / (2*x_std^2));   
    gtmp = gtmp / sum(gtmp);
    g = g * (max(gtmp) / max(g));
    plot(xg,g,'r','linewidth',2)
    title(sprintf('Mean = %0.3f  Std = %0.3f', x_mean, x_std),'FontSize',14);
    xlabel('x error (pixels)','FontSize',14)
    ylabel('Probability','FontSize',14)
    set(gca,'FontSize',14);
    set(gca,'XLim',[-alim alim])
    grid on
    set(gcf,'PaperPositionMode','auto')
    %svfile = [pathdir sprintf('/reproj_xhist_all_homog_poly%0.2d.png',porder)];
    %print('-dpng',svfile)

    figure(11);
    cla;
    hold on
    bar(ybins, hy/sum(hy));
    y_mean = mean(ey);
    y_std = std(ey);
    yg = -alim:0.1:alim;
    g = exp(-(yg-y_mean).^2 / (2*y_std^2));
    gtmp = exp(-(ybins-y_mean).^2 / (2*y_std^2));   gtmp = gtmp / sum(gtmp);
    g = g * (max(gtmp) / max(g));
    plot(yg,g,'r','linewidth',2)
    title(sprintf('Mean = %0.3f  Std = %0.3f', y_mean, y_std),'FontSize',14);
    xlabel('y error (pixels)','FontSize',14)
    ylabel('Probability','FontSize',14)
    set(gca,'FontSize',14);
    set(gca,'XLim',[-alim alim])
    grid on
    % set(gcf,'PaperPositionMode','auto')
    % svfile = [pathdir sprintf('/reproj_yhist_all_homog_poly%0.2d.png',porder)];
    % print('-dpng',svfile)

    
    
    % Plot error and function of distance from principal point
    %===========================================================
    figure(14)
    etotal = sqrt(ex.^2 + ey.^2);
    rad = rad(:);
    rlim = 50 * ceil(max(rad)/50);
    cla;
    hold on
    plot(rad, etotal,'b.')

    % Run a 20 frame sliding window average, each 5 pixels
    rbins = 0:10:max(rad);
    err_med_r = -1 * ones(1,length(rbins));
    err_mean_r = -1 * ones(1, length(rbins));
    for i = 1:length(rbins)
        rstart = rbins(i)-25;
        rend = rbins(i)+25;
        ind = find(rad >= rstart & rad <= rend);
        if ~isempty(ind)
            vals = etotal(ind);
            err_med_r(i) = median(vals);
            err_mean_r(i) = mean(vals);
        end
    end
    ind = find(err_med_r >= 0);
    rbins = rbins(ind);
    err_med_r = err_med_r(ind);
    err_mean_r = err_mean_r(ind);
    plot(rbins, err_med_r, '-r','linewidth',2)
    plot(rbins, err_mean_r, '-k','linewidth',2)
    set(gca,'XLim',[0 rlim],'XTick',0:50:rlim)

    ylabel('Total error (pixels)','FontSize',14)
    xlabel('Radius from Principal Point (pixels)','FontSize',14)
    legend('Obs','Median','Mean','Location','NorthWest')
    set(gca,'FontSize',14);
    grid on
    % set(gcf,'PaperPositionMode','auto')
    % svfile = [pathdir sprintf('/reproj_vsrad_homog_poly%0.2d.png',porder)];
    % print('-dpng',svfile)
    
    
    
    
end



% 
%     
% % Plot all the errors together (points)
% figure(8)
% cla;
% hold on
% plot(err_u(:), err_v(:), 'bx','markersize',5)
% plot(mean(err_u(:)),mean(err_v(:)),'rs','linewidth',3)
% phi = linspace(0,2*pi,100);
% plot(0.5*cos(phi),0.5*sin(phi),'k')
% plot(cos(phi),sin(phi),'k')
% axis equal
% alim = ceil(2*max(abs([err_u(:); err_v(:)])))/2;
% set(gca,'XLim',[-alim alim],'YLim',[-alim alim]);
% set(gca,'XTick',-alim:0.5:alim,'YTick',-alim:0.5:alim);
% set(gca,'FontSize',14)
% title(sprintf('Mean = (%0.2f , %0.2f)', mean(err_u(:)), mean(err_v(:))),'FontSize',14)
% xlabel('x error (pixels)','FontSize',14)
% ylabel('y error (pixels)','FontSize',14)
% set(gcf,'PaperPositionMode','Auto')
% grid on
% box on    
% % set(gcf,'PaperPositionMode','auto')
% % svfile = [pathdir sprintf('/reproj_xy_all_homog_poly%0.2d.png',porder)];
% % print('-dpng',svfile)



% 
% 
% figure(12);
% cla
% etotal = sqrt(ex.^2 + ey.^2);
% alim = ceil(4 * max(etotal)) / 4;
% xbins = (0:0.125:alim) + 0.125/2;
% h = hist(etotal, xbins);
% bar(xbins, h/length(etotal));
% xlabel('Total error (pixels)','FontSize',14)
% ylabel('Probability','FontSize',14)
% set(gca,'FontSize',14);
% set(gca,'XLim',[0 alim])
% grid on
% % set(gcf,'PaperPositionMode','auto')
% % svfile = [pathdir sprintf('/reproj_hist_all_homog_poly%0.2d.png',porder)];
% % print('-dpng',svfile)
% 
% 
% 
% 
% 
% % Plot the reprojection error vs. distance from the principal point
% %=========================================================================
% figure(14)
% rlim = 50 * ceil(max(rad)/50);
% cla;
% hold on
% plot(rad, etotal,'b.')
% 
% % Run a 20 frame sliding window average, each 5 pixels
% rbins = 0:10:max(rad);
% err_med_r = -1 * ones(1,length(rbins));
% err_mean_r = -1 * ones(1, length(rbins));
% for i = 1:length(rbins)
%     rstart = rbins(i)-25;
%     rend = rbins(i)+25;
%     ind = find(rad >= rstart & rad <= rend);
%     if ~isempty(ind)
%         vals = etotal(ind);
%         err_med_r(i) = median(vals);
%         err_mean_r(i) = mean(vals);
%     end
% end
% ind = find(err_med_r >= 0);
% rbins = rbins(ind);
% err_med_r = err_med_r(ind);
% err_mean_r = err_mean_r(ind);
% plot(rbins, err_med_r, '-r','linewidth',2)
% plot(rbins, err_mean_r, '-k','linewidth',2)
% set(gca,'XLim',[0 rlim],'XTick',0:50:rlim)
% 
% ylabel('Total error (pixels)','FontSize',14)
% xlabel('Radius from Principal Point (pixels)','FontSize',14)
% legend('Obs','Median','Mean','Location','NorthWest')
% set(gca,'FontSize',14);
% grid on
% % set(gcf,'PaperPositionMode','auto')
% % svfile = [pathdir sprintf('/reproj_vsrad_homog_poly%0.2d.png',porder)];
% % print('-dpng',svfile)
% 

