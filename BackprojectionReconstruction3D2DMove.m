% y axis - horizontal printer move
% x axis - vetical printer move
% z axis - towards target

% Backprojection based off of grid center (assume NxN square grid)
%close all; clc; clear;
disp('Reconstruction Starting'); tic;

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'YZ_2D_20250728_143347.mat';
close all; disp('Reconstruction Starting'); tic;

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'YZ_2D_20250730_171530.mat';
load(dataFile)   % contains variable `recs`
[~, baseName] = fileparts(dataFile);

vtrigU_ants_location;
[Xgrid,Ygrid,Zgrid]=meshgrid(xgrid,ygrid,zgrid);

src = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);
src2 = permute(src,[3,2,4,1]);
c = physconst('lightspeed'); %(m/s)
N_freq = length(freq);
Nfft = 2^(ceil(log2(size(freq,2)))+1);

%% dynamic grid collapse (collapse the smallest axis) ----
% lengths of each axis
lens = [numel(xgrid), numel(ygrid), numel(zgrid)];  
% pick the axis with the fewest samples
[~, sliceDim] = min(lens);            

% collapse that axis to its first value
gv = {xgrid, ygrid, zgrid};
gv{sliceDim} = gv{sliceDim}(1);       
[xgrid, ygrid, zgrid] = deal(gv{:});

% determine which two axes remain for rows/cols
denseDims = setdiff(1:3, sliceDim);  
gridA     = gv{denseDims(1)};  % rows
gridB     = gv{denseDims(2)};  % cols
Ny = numel(gridA);
Nx = numel(gridB);

%% pre‑allocate reconstruction stack -------------
nRecs   = size(recs,3);                 % how many sweeps in .mat
y_accum = zeros(Ny, Nx, nRecs, 'single');   % or 'double'
y_cart_sum = zeros(Ny, Nx, 'single');
toDB = @(M) 20*log10(abs(M)+eps) - max(20*log10(abs(M(:))+eps)); %db helper function


%% Back Projection Loop --------------------------
for i=1:nRecs
X = recs(:,:,i);
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
padsig = [padsig((lnconv-1)/2:-1:1),padsig,padsig(end:-1:end-(lnconv-1)/2+1)]; 
padsig = conv(padsig,c2,'valid');        
f_res = padsig>thresh;

%Remove resonant frequencies
X = X .* (1-f_res);  
%convert to complex time domain signal
x = ifft(X,Nfft,2);
       
y_cart = reshape(H2*reshape(X,[],1),size(Xgrid));

% pull out the slice (will be linear if one dim was collapsed)
temp_y = squeeze(y_cart);    % → [Ny×Nx] matrix every time

% now safe to accumulate
y_accum(:,:,i) = temp_y;

y_cart_sum = y_cart_sum + temp_y; % accumulate for reconstruction

a = sprintf('Image %d/%d Proccessed',i, numel(recs(1,1,:)) ); disp(a);

end %% End of Backprojection
%% decide which plane we have --------------------
axisLabels = 'XYZ';            % 1:X  2:Y  3:Z

if     numel(xgrid) <= 10
        sliceName = 'YZ';  rows = ygrid; cols = zgrid;
elseif numel(ygrid) <= 10
        sliceName = 'XZ';  rows = xgrid; cols = zgrid;
else
        sliceName = 'XY';  rows = xgrid; cols = ygrid;

        sliceName = 'YX';  rows = xgrid; cols = ygrid;
end

%% Improved Image -------------------------------------------
C_imp = y_cart_sum;     % or rssq(y_accum,3) for RSS average
Norm_C_imp = C_imp ./ max(abs(C_imp(:)));
plotSlice(Norm_C_imp, rows, cols, sliceName, 'Improved', baseName)

%% Unimproved (first record) ---------------------------------
C_raw = y_accum(:,:,ceil(nRecs/2) ); %return grid center image
Norm_C_raw = C_raw ./ max(abs(C_raw(:)));
plotSlice(Norm_C_raw, rows, cols, sliceName, 'Unimproved', baseName)

%% Helper Functions
% function C = getSlice(Y, sliceDim, Ny, Nx)
%     % Build an index that fixes sliceDim→1, leaves the other dims free
%     idx = repmat({':'}, 1, ndims(Y));
%     idx{sliceDim} = 1;
% 
%     % Extract the slab—this may come out as 1×Ny×Nx, Ny×1×Nx, Ny×Nx×1, or even 1×N-vector
%     slab = Y(idx{:});
% 
%     % Flatten to a column vector
%     vec = slab(:);
% 
%     % Sanity check element count
%     if numel(vec) ~= Ny * Nx
%         error('getSlice: element mismatch — got %d, expected %d', ...
%                numel(vec), Ny * Nx);
%     end
% 
%     % Reshape into Ny-by-Nx
%     C = reshape(vec, [Ny, Nx]);
% end



function plotSlice(C, rows, cols, sliceName, tag, baseName)
    % ensure we’re working with a 2‑D slice
    C2 = squeeze(C);                        
    assert(ndims(C2)==2, 'plotSlice: expected 2‑D input, got %d‑D', ndims(C2));

    % convert to dB
    Cdb = 20*log10(abs(C2) + eps);

    figure
    % swap rows/cols via permute instead of .' on an N‑D
    imagesc(cols, rows, permute(Cdb, [2,1]))
    axis xy equal tight
    colormap turbo; colorbar

    % (the rest stays the same…)
    lim = caxis; lim(1)=floor(lim(1)/5)*5; lim(2)=ceil(lim(2)/5)*5;
    caxis(lim)
    xlabel(sprintf('%c [m]', sliceName(2)))
    ylabel(sprintf('%c [m]', sliceName(1)))
    title(sprintf('%s Power Slice (%s)', sliceName, tag))

    outFile = sprintf('%s_%s_%s.png', baseName, sliceName, tag);
    exportgraphics(gcf, outFile, 'Resolution', 300)
end
toc
