
%load reconstructed file
load('reconstructed_ball.mat');

figure();
nRecs = numel(recs(1,1,:));

for i=1:nRecs:1
%% Plot Y-Z Slice
        if and(min([length(ygrid),length(zgrid)])>2,length(xgrid)<=10)
            subplot(nRecs/5,nRecs/5, i)
            y_yz = 20*log10(rssq(y_cart(:,find(xgrid>=xgrid(1),1):find(xgrid>=xgrid(end),1),:),2));
            ax=pcolor(squeeze(Ygrid(:,1,:)),squeeze(Zgrid(:,1,:)),squeeze(y_yz));
            set(ax,'EdgeColor', 'none');
            xlabel('Y [m]'); ylabel('Z [m]');
            title(sprintf("YZ Scan %d", i));
            filename = sprintf('YZ_2D_%s', dateTag);

        end
%% Plot X-Z Slice
        if and(min([length(xgrid),length(zgrid)])>2,length(ygrid)<=10)
            y_xz = 20*log10(rssq(y_cart(find(ygrid>=ygrid(1),1):find(ygrid>=ygrid(end),1),:,:),1));
            figure(fig(3));ax=pcolor(squeeze(Xgrid(1,:,:)),squeeze(Zgrid(1,:,:)),squeeze(y_xz));
            set(ax,'EdgeColor', 'none');
            xlabel('X [m]'); ylabel('Z [m]');
            filename = sprintf('XZ_2D_%s', dateTag);
            % if first_iter 
            %     set(gca,'NextPlot','replacechildren');
            %     title('xz view');xlabel('x');ylabel('z');daspect([1,1,1]);%caxis([-20,20]);
            % end
        end

%% Plot X-Y Slice
        if and(min([length(xgrid),length(ygrid)])>2,length(zgrid)<=10)
           % y_xy = 20*log10(rssq(y_cart(:,:,find(zgrid>=zgrid(1),1):find(zgrid>=zgrid(end),1)),2));
            figure(fig(2));ax = pcolor( squeeze(Xgrid(:,:,1)), squeeze(Ygrid(:,:,1)), ...
             20*log10(rssq(y_cart(:,:,find(zgrid>=zgrid(1),1):find(zgrid>=zgrid(end),1)),3)) );
            filename = sprintf('XY_2D_%s', dateTag);
            % if first_iter
            %     set(gca,'NextPlot','replacechildren');
            %     title('xy view');xlabel('x');ylabel('y');daspect([1,1,1]);%caxis([-20,20]);
            % end
        end
end %end for loop
