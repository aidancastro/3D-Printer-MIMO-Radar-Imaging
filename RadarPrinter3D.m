%add DOTNET assembly
clc;
close all;
clear;
NET.addAssembly([getenv('programfiles'),'\Vayyar\vtrigU\bin\vtrigU.CSharp.dll']);
%load('Calibration_YZ_2D_20250730_162621.mat'); tic;

%calibration_X = recs;

%export file info
extension = '.mat'; export_path = "C:\Users\13523\Desktop\Vakalis research\Figure Holding Space";
dateTag   = datestr(datetime('now'),'yyyymmdd_HHMMss'); 

%% Dynamic Square Snake scanning grid 
NxN=5; xstep=40; zstep=40; gridCenter=[100,100]; 
printer_positions=zeros(NxN^2,2);printer_offsets = zeros(NxN^2,2); ctr=ceil(NxN/2); idx=1;
for row=1:NxN
    z_off=(row-ctr)*zstep;
    cols = (mod(row,2)==1)*(1:NxN) + (mod(row,2)==0)*(NxN:-1:1);
    for col=cols
        x_off=(col-ctr)*xstep;
        printer_positions(idx,:) = gridCenter+[x_off z_off];
        printer_offsets(idx,:)   = [x_off z_off];   
        idx=idx+1;
    end
end
printer_offsets = printer_offsets .* 0.001;
%% Radar Settings

%init the device
vtrigU.vtrigU.Init();
%Set setting structure
mySettings.freqRange.freqStartMHz = 62*1000;
mySettings.freqRange.freqStopMHz = 69*1000;
mySettings.freqRange.numFreqPoints = 30;
mySettings.rbw_khz = 40;
mySettings.txMode = vtrigU.TxMode.LOW_RATE;
vtrigU.vtrigU.ApplySettings(mySettings.freqRange.freqStartMHz,mySettings.freqRange.freqStopMHz, mySettings.freqRange.numFreqPoints, mySettings.rbw_khz, mySettings.txMode );

%Validate settings
mySettings = vtrigU.vtrigU.GetSettings();
vtrigU.vtrigU.ValidateSettings(mySettings);

%% Get TxRx Pairs and convert to complex matrix
ants_pairs = vtrigU.vtrigU.GetAntennaPairs(mySettings.txMode);
TxRxPairs = zeros(ants_pairs.Length,2);
for ii = 1: ants_pairs.Length
    TxRxPairs(ii,1) = double(ants_pairs(ii).tx);
    TxRxPairs(ii,2) = double(ants_pairs(ii).rx);
end

%Get used frequencies in Hz
freq = double(vtrigU.vtrigU.GetFreqVector_MHz())*1e6;

%get antennas locations from script
vtrigU_ants_location;

%Define constants
N_txrx = size(TxRxPairs,1);
N_freq = length(freq);
c = 3e8; %SOL (m/s)

Nfft = 2^(ceil(log2(size(freq,2)))+1);
Ts = 1/Nfft/(freq(2)-freq(1)+1e-16); %Avoid nan checks
time_vec = (0:Ts:Ts*(Nfft-1));
dist_vec = time_vec*1.5e8; %distance in meters

%% Voxels
Nx = 30; Ny = 30; Nz = 20;
xgrid = linspace(-0.1,0.4, Nx);
zgrid = linspace(0.1,0.4,Nz);
ygrid = linspace(-0.1, 0.4, Ny);
[Xgrid,Ygrid,Zgrid]=meshgrid(xgrid,ygrid,zgrid);

src = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);
src2 = permute(src,[3,2,4,1]);

% Use Radar equation to find total loss
% Pe/Ps = (Gtx*Grx*RCS*lambda^2)/((4pi)^3*|Rtx|^2*|Rrx|^2)
c = physconst('lightspeed'); %(m/s)

%% Beam Steering
Rvec = src2- VtrigU_ants_location;
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

%% Start Recording

nfig = 4; fig = zeros(1,nfig);
for ii = 1:nfig
    fig(ii) = figure();
end

threeDPrinter("open");

threeDPrinter("move", printer_positions(1,1), printer_positions(1,2)); pause(10);
prev_pos=[0 0];

for kk = 1:length(printer_positions)
        pos=printer_positions(kk,:);
    if pos(1)==prev_pos(1)
        threeDPrinter("move",pos(1),pos(2)); pause(9); %vertical movement
    else
        threeDPrinter("move",pos(1),pos(2)); pause(2); %horizontal movement
    end
    prev_pos=pos;

        vtrigU.vtrigU.Record();
        rec = double(vtrigU.vtrigU.GetRecordingResult(vtrigU.SignalCalibration.DEFAULT_CALIBRATION));

        X = rec;
       
        %Convert and reshape result to matlab complex matrix
        smat_size = [size(TxRxPairs,1),size(freq,2),2];
        my_perms = [1,3,2];
        X = reshape(X,smat_size(my_perms));
        X = ipermute(X,my_perms);
        X = X(:,:,1)+ 1j*X(:,:,2);
        recs(:,:,kk) = X;
        
            if kk == 1
            
            figure(fig(1));
            clf; hold on;
            for ll = 1:smat_size(1)
                plot(freq,20*log10(abs(X(ll,:))));
            end
            plot(freq,20*log10(rssq(X,1))-10*log10(smat_size(1)),'k.','LineWidth',4);
            title('Channel response (per channel and rms average)')
            
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
            end  
        
        %Remove resonant frequencies
        X = X .* (1-f_res);
        %X = X - calibration_X(:,:,kk); %calibration
        
        %convert to complex time domain signal
        x = ifft(X,Nfft,2);
       
        y_cart = reshape(H2*reshape(X,[],1),size(Xgrid));

        %Create and show power delay profile - non coherent summation
        PDP = rssq(x,1);
        figure(fig(3));plot(dist_vec,20*log10(abs(PDP./max(abs(PDP)))));
        
%% Scatter Plot 
        if min([length(xgrid),length(ygrid),length(zgrid)])>2
            th = abs(y_cart)>1;
            figure(fig(2));
            scatter3(Xgrid(th),Ygrid(th),Zgrid(th),20*log10(abs(y_cart(th))),20*log10(abs(y_cart(th))))
            filename = sprintf('YZ_2D_%s', dateTag);

        end
%% Plot Y-Z Slice
        if and(min([length(ygrid),length(zgrid)])>2,length(xgrid)<=10)
            y_yz = 20*log10(rssq(y_cart(:,find(xgrid>=xgrid(1),1):find(xgrid>=xgrid(end),1),:),2));
            figure(fig(2));ax=pcolor(squeeze(Ygrid(:,1,:)),squeeze(Zgrid(:,1,:)),squeeze(y_yz));
            set(ax,'EdgeColor', 'none');
            xlabel('Y [m]'); ylabel('Z [m]');
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
        %first_iter=false;
        pause(3)
end %% end scan loop

save(filename, 'y_cart', 'xgrid','ygrid','zgrid','freq','TxRxPairs','NxN','xstep','zstep','gridCenter','printer_offsets');

threeDPrinter("home");
threeDPrinter("close");

toc;

