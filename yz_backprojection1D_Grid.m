% y axis - horizontal printer move
% x axis - vetical printer move
% z axis - towards target
%%
clc; clear;

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'radar_scan_dataXY1.mat';
load(dataFile)   % contains variable `recs`
[~, baseName] = fileparts(dataFile);

vtrigU_ants_location;
[Xgrid,Ygrid,Zgrid]=meshgrid(xgrid,ygrid,zgrid);

src = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);
src2 = permute(src,[3,2,4,1]);

c = physconst('lightspeed'); %(m/s)
N_freq = length(freq);
Nfft = 2^(ceil(log2(size(freq,2)))+1);
y_cart_sum = zeros(size(xgrid));
single_scan = X(:,:,1); %first image for reference

%% Back Projection Loop 
for i=1:10
X = recs(:,:,i);

Rvec = src2-(VtrigU_ants_location + [0 0.02*(i-1) 0]);
Rmag = rssq(Rvec,2);
Rtheta = atan2(rssq(Rvec(:,1:2,:,:),2),Rvec(:,3,:,:));
Rphi = atan2(Rvec(:,2,:,:),Rvec(:,1,:,:));
Sphase = 2*pi*Rmag.*freq/c; %Electrical Length in Radians
RCS = 1; %m^2
lambda = c./freq; csf = sqrt(RCS).*lambda./((4*pi).^(3/2));
Smag = 10^(5.8/20)*RadiationPattern(Rtheta,Rphi)./Rmag;

H2 = zeros(length(TxRxPairs),length(freq),1,length(src2));
for ii = 1:length(TxRxPairs)
    tx = TxRxPairs(ii,1); rx = TxRxPairs(ii,2);
    H2(ii,:,:,:) = 1./(csf.*Smag(tx,:,:,:).*Smag(rx,:,:,:).*...
              exp(-1j.*(Sphase(tx,:,:,:)+Sphase(rx,:,:,:))));
end
H2 = reshape(permute(H2,[4,1,2,3]),length(src2),[]); %xyz x txrx x freq

 % Identify resonant frequencies 
thresh = 3;
lnconv = min(max(floor(N_freq/8)*2+1,floor(50/(freq(2)-freq(1)))*2+1),...
         floor(3*N_freq/8)*2+1); %conv length between 1/4 and 3/4 N_freq
c2 = -ones(lnconv,1)/(lnconv-1);
c2((lnconv+1)/2) = 1;
padsig = 20*log10(rssq(X,1));
padsig = [padsig((lnconv-1)/2:-1:1),padsig,padsig(end:-1:end-(lnconv-1)/2+1)]; 
padsig = conv(padsig,c2,'valid');        
f_res = padsig>thresh;

%Remove resonant frequencies
X = X .* (1-f_res);  
%convert to complex time domain signal
x = ifft(X,Nfft,2);
       
y_cart = reshape(H2*reshape(X,[],1),size(Xgrid));

y_cart_sum = y_cart_sum + y_cart;

disp(i)
end %end for 

<<<<<<< Updated upstream
<<<<<<< Updated upstream
 %% Plot Y-Z Slice
        if and(min([length(ygrid),length(zgrid)])>2,length(xgrid)<=2)
            y_yz = 20*log10(rssq(y_cart_sum(:,find(xgrid>=xgrid(1),1):find(xgrid>=xgrid(end),1),:),2));
            figure();ax=pcolor(squeeze(Ygrid(:,1,:)),squeeze(Zgrid(:,1,:)),squeeze(y_yz));
            set(ax,'EdgeColor', 'none');
            % if first_iter 
            %     set(gca,'NextPlot','replacechildren');
            %     title('yz view');xlabel('y');ylabel('z');daspect([1,1,1]);%caxis([-20,20]); 
            % end
            outname = sprintf('%s_yz.png',baseName);
            exportgraphics(gcf, outName, 'Resolution', 300);
        end
        %% Plot X-Z Slice
        if and(min([length(xgrid),length(zgrid)])>2,length(ygrid)<=2)
            y_xz = 20*log10(rssq(y_cart_sum(find(ygrid>=ygrid(1),1):find(ygrid>=ygrid(end),1),:,:),1));
            figure();ax=pcolor(squeeze(Xgrid(1,:,:)),squeeze(Zgrid(1,:,:)),squeeze(y_xz));
            set(ax,'EdgeColor', 'none');
            % if first_iter 
            %     set(gca,'NextPlot','replacechildren');
            %     title('xz view');xlabel('x');ylabel('z');daspect([1,1,1]);%caxis([-20,20]);
            % end
        end
=======
=======
>>>>>>> Stashed changes
%% Plot YZ-Slice
 if  min([numel(ygrid) numel(zgrid)]) > 2 && numel(xgrid) <= 2

    % 1. Collapse along X and convert to dB
    y_yz = 20*log10( rssq(y_cart_sum, 2) );   % size = Ny × Nz
    y_yz = y_yz - max(y_yz(:)); % peak-normalize
    % 2. Plot
    figure;
    p = pcolor( squeeze(Ygrid(:,1,:)), ...
                squeeze(Zgrid(:,1,:)), ...
                squeeze(y_yz) );
    set(p,'EdgeColor','none');
    shading interp;             % smoother look
    colormap(turbo);            % colour‑blind‑safe

    % 3. Dynamic colour limits: top 40 dB
    peak = max(y_yz(:));
    caxis([-40 0]);

    cb = colorbar;  ylabel(cb,'|χ| [dB]');

    axis equal tight;           % square pixels, no margins
    xlabel('Y [m]');  ylabel('Z [m]');
    title('Y–Z Power Slice');
 end
 %% Plot XZ-Slice
 if  min([numel(xgrid) numel(zgrid)]) > 2 && numel(ygrid) <= 2
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

    % 1. Collapse along X and convert to dB
    y_xz = 20*log10( rssq(y_cart_sum, 2) );   % size = Ny × Nz
    y_xz = y_xz - max(y_xz(:)); % peak-normalize
    % 2. Plot
    figure;
    p = pcolor( squeeze(Xgrid(:,1,:)), ...
                squeeze(Zgrid(:,1,:)), ...
                squeeze(y_xz) );
    set(p,'EdgeColor','none');
    shading interp;             % smoother look
    colormap(turbo);            % colour‑blind‑safe

    % 3. Dynamic colour limits: top 40 dB
    peak = max(y_xz(:));
    caxis([-40 0]);

    cb = colorbar;  ylabel(cb,'|χ| [dB]');

    axis equal tight;           % square pixels, no margins
    xlabel('X [m]');  ylabel('Z [m]');
    title('X–Z Power Slice');
 end
 %% Plot XY-Slice
 if  min([numel(xgrid) numel(ygrid)]) > 2 && numel(zgrid) <= 2

    % 1. Collapse along X and convert to dB
    y_xy = 20*log10( rssq(y_cart_sum, 2) );   % size = Ny × Nz
    y_xy = y_xy - max(y_xy(:)); % peak-normalize
    % 2. Plot
    figure;
    p = pcolor( squeeze(Xgrid(:,1,:)), ...
                squeeze(Ygrid(:,1,:)), ...
                squeeze(y_xy) );
    set(p,'EdgeColor','none');
    shading interp;             % smoother look
    colormap(turbo);            % colour‑blind‑safe

    % 3. Dynamic colour limits: top 40 dB
    peak = max(y_xy(:));
    caxis([-40 0]);

    cb = colorbar;  ylabel(cb,'|χ| [dB]');

    axis equal tight;           % square pixels, no margins
    xlabel('X [m]');  ylabel('Y [m]');
    title('X-Y Power Slice');
 end