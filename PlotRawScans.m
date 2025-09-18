%script is used to plot raw recordings using the same steering matrix
%expects inputs recs, xgrid,ygrid,zgrid,Xgrid,Ygrid,Zgrid,
%vTrigU_antenna_locations,

vtrigU_ants_location;

dataFile = 'ball_Scan_20250917_112419.mat';
load(dataFile)   % recs, freq, xgrid, ygrid, zgrid, printer_offsets
[~, baseName] = fileparts(dataFile);

src = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);
src2 = permute(src,[3,2,4,1]);

nFreq = numel(recs(1,:,1));
freq = linspace(62e9,69e9, nFreq);

% Use Radar equation to find total loss
% Pe/Ps = (Gtx*Grx*RCS*lambda^2)/((4pi)^3*|Rtx|^2*|Rrx|^2)
c = physconst('lightspeed'); %(m/s)

%% Beam Steering
Rvec = src2-VtrigU_ants_location;
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
size(H2)

y_cart = reshape(H2.*reshape(recs,[],1),size(Xgrid)); %reshape for display while steering

%%start loop
for i=1:nRecs
figure(i);

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
