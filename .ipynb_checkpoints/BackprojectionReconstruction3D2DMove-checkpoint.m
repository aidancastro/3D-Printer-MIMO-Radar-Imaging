% y axis - horizontal printer move
% x axis - vetical printer move
% z axis - towards target

% Backprojection based off of grid center (assume NxN square grid)
%close all; clc; clear;

dataFile = '3D_20250730_193754.mat';
close all; disp('Reconstruction Starting'); tic;

%% ---------------- USER PARAMS ----------------------------------------
load(dataFile)   % contains variable `recs`
[~, baseName] = fileparts(dataFile);

vtrigU_ants_location;
[Xgrid,Ygrid,Zgrid]=meshgrid(xgrid,ygrid,zgrid);

src = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);
src2 = permute(src,[3,2,4,1]);
c = physconst('lightspeed'); %(m/s)
N_freq = length(freq);
Nfft = 2^(ceil(log2(size(freq,2)))+1);

%% pre‑allocate reconstruction stack -------------

%length of each voxel axis
Ny = numel(xgrid);
Nx = numel(ygrid);
Nz = numel(zgrid);

nRecs   = size(recs,3);                 % how many sweeps in .mat
y_accum = zeros(Ny, Nx, Nz, nRecs, 'single');   % or 'double'
y_cart_sum = zeros(Ny, Nx, Nz,'single');
toDB = @(M) 20*log10(abs(M)+eps) - max(20*log10(abs(M(:))+eps)); %db helper function


%% Back Projection Loop --------------------------
for i=1:nRecs
X = squeeze(recs(:,:,i));
Rvec = src2-(VtrigU_ants_location + [printer_offsets(i,1) printer_offsets(i,2) 0]); %printer_offsets..  y(R->L) x(up) from vtrigU
Rmag = rssq(Rvec,2);
Rtheta = atan2(rssq(Rvec(:,1:2,:,:),2),Rvec(:,3,:,:));
Rphi = atan2(Rvec(:,2,:,:),Rvec(:,1,:,:));
Sphase = 2*pi*Rmag.*freq/c; %Electrical Length in Radians
RCS = 1; %m^2
lambda = c./freq; csf = sqrt(RCS).*lambda./((4*pi).^(3/2));
Smag = 10^(5.8/20)*RadiationPattern(Rtheta,Rphi)./Rmag;

% beam steering matrix
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
%padsig = padsig(:).';  % ensure row vector

padsig = [padsig((lnconv-1)/2:-1:1),padsig,padsig(end:-1:end-(lnconv-1)/2+1)]; 
padsig = conv(padsig,c2,'valid');        
f_res = padsig>thresh;

%Remove resonant frequencies
X = X .* (1 - f_res);
%convert to complex time domain signal
x = ifft(X,Nfft,2);
       
y_cart = reshape(H2*reshape(X,[],1), numel(xgrid), numel(ygrid), numel(zgrid));
temp_y = y_cart;

% accumulate
y_accum(:,:,:,i) = temp_y;

y_cart_sum = y_cart_sum + temp_y; % accumulate for reconstruction

a = sprintf('Image %d/%d Proccessed',i, numel(recs(1,1,:)) ); disp(a);

end %% End of Backprojection

%% Plot Images -------------------------------------------
C_imp = y_cart_sum;     % improved from nRecs images
C_raw = y_accum(:,:,:,ceil(nRecs/2) ); %return grid center image

plotRadar3DCompare(C_imp,C_raw, xgrid, ygrid, zgrid, 'isosurface');

%% Helper Functions
function plotRadar3DCompare(CubeImproved, CubeUnimproved, xgrid, ygrid, zgrid, style)
% plotRadar3DCompare  Plots Improved vs Unimproved 3D radar reconstructions
%
% Inputs:
%   CubeImproved, CubeUnimproved - 3D voxel cubes (Ny x Nx x Nz)
%   xgrid, ygrid, zgrid - coordinate vectors
%   style - 'isosurface', 'slice', 'scatter', 'surf'

[X, Y, Z] = meshgrid(xgrid, ygrid, zgrid);

% Normalize magnitudes
CubeImpNorm  = abs(CubeImproved)  ./ max(abs(CubeImproved(:)));
CubeUnimpNorm = abs(CubeUnimproved) ./ max(abs(CubeUnimproved(:)));

figure;
for k = 1:2
    subplot(1,2,k)
    if k == 1
        CubeNorm = CubeUnimpNorm;
        titleStr = 'Unimproved';
    else
        CubeNorm = CubeImpNorm;
        titleStr = 'Improved';
    end
    
    switch lower(style)
        case 'isosurface'
            thresh = 0.5; % ~ -6 dB
            p = patch(isosurface(X, Y, Z, CubeNorm, thresh));
            isonormals(X, Y, Z, CubeNorm, p);
            set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
            camlight; lighting gouraud;
            
        case 'slice'
            xs = xgrid(round(end/2));
            ys = ygrid(round(end/2));
            zs = zgrid(round(end/2));
            slice(X, Y, Z, CubeNorm, xs, ys, zs);
            shading interp;
            
        case 'scatter'
            mask = CubeNorm > 0.5;
            scatter3(X(mask), Y(mask), Z(mask), 10, CubeNorm(mask), 'filled');
            
        case 'surf'
            idx = round(length(zgrid)/2);
            surf(X(:,:,idx), Y(:,:,idx), CubeNorm(:,:,idx));
            shading interp;
            
        otherwise
            error('Unknown style. Use isosurface, slice, scatter, or surf');
    end
    
    title([titleStr ' Reconstruction']);
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    axis tight; grid on; view(3);
end



end

toc
