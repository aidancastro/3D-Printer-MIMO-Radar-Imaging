
close all;
%load reconstructed file
%expects x(X)grid,y(Y)grid,z(Z)grid, y_accum, 
load('reconstructed_ball_Scan_20250917_111655_YX.mat');

%figure('Position',[100 100 1400 900]); % width x height in pixels
nRecs = numel(y_accum(1,1,:));

for i=1:nRecs
figure(i);
y_cart = y_accum(:,:,i);
%% Plot Y-Z Slice
        if and(min([length(ygrid),length(zgrid)])>2,length(xgrid)<=10)
            %subplot(3,1, j);
            y_yz =20*log10(rssq(y_cart,3));
            ax=pcolor(squeeze(Ygrid),squeeze(Zgrid),squeeze(y_yz));
            set(ax,'EdgeColor', 'none');
            xlabel('Y [m]'); ylabel('Z [m]');
            title(sprintf("YZ Scan %d", i));
            

        end
%% Plot X-Z Slice
        if and(min([length(xgrid),length(zgrid)])>2,length(ygrid)<=10)
            %subplot(5,1, j);
            y_xz = 20*log10(rssq(y_cart,3));
            ax=pcolor(squeeze(xgrid(1,:,:)),squeeze(zgrid(1,:,:)),squeeze(y_xz));
            set(ax,'EdgeColor', 'none');
            xlabel('X [m]'); ylabel('Z [m]');
            
            % if first_iter 
            %     set(gca,'NextPlot','replacechildren');
            %     title('xz view');xlabel('x');ylabel('z');daspect([1,1,1]);%caxis([-20,20]);
            % end
        end

%% Plot X-Y Slice
        if and(min([length(xgrid),length(ygrid)])>2,length(zgrid)<=10)
            %subplot(5,1, j);
            y_xy = 20*log10(rssq(y_cart,3));
            ax = pcolor( squeeze(xgrid(:,:,1)), squeeze(ygrid(:,:,1)), ...
             20*log10(rssq(y_cart(:,:,find(zgrid>=zgrid(1),1):find(zgrid>=zgrid(end),1)),3)) );
            
            % if first_iter
            %     set(gca,'NextPlot','replacechildren');
            %     title('xy view');xlabel('x');ylabel('y');daspect([1,1,1]);%caxis([-20,20]);
            % end
        end
        fprintf("completed plot %d", i);
end %end for loop

