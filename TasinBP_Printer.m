clear all
close all;

%% initial measurements from the IMAGEVK Radar
load("Gun_Scan_20250908_201059.mat");
nRecs = size(recs,3);                      % safer
c = 299792458;

% antenna locations
load("AntennaLocations.mat");              % TX_x, TX_y, RX_x, RX_y expected

% basic sizes
nx  = numel(xgrid);
ny  = numel(ygrid);
nTx = numel(TX_x);
nRx = numel(RX_y);
nF  = numel(freq);

V     = nx*ny;
Nobs  = nTx*nRx*nF;
Im    = complex(zeros(V,1));               % keep vector then reshape once

% z offset (down-range)
z_offset = single(0.5);

% voxel grid once (nx x ny)
[Xg, Yg] = ndgrid(single(xgrid), single(ygrid));   % Xg,Yg: nx x ny

% Pre-allocate backprojected antenna locations
BP_TX_x = ones(nRecs, nTx, 'single');
BP_TX_y = ones(nRecs, nTx, 'single');
BP_RX_x = ones(nRecs, nRx, 'single');
BP_RX_y = ones(nRecs, nRx, 'single');

% Build backprojected antenna locations per scan (printer_offsets = [yoff xoff])
for i = 1:nRecs
    BP_TX_x(i,:) = single(TX_x(:)).' - single(printer_offsets(i,2));
    BP_TX_y(i,:) = single(TX_y(:)).' - single(printer_offsets(i,1));
    BP_RX_x(i,:) = single(RX_x(:)).' - single(printer_offsets(i,2));
    BP_RX_y(i,:) = single(RX_y(:)).' - single(printer_offsets(i,1));
end

%% Main accumulation over scans
for i = 1:nRecs
    recording = recs(:,:,i);               % (nTx*nRx, nF) or (nPairs, nF)
    % If you truly have TX×RX flattened in TX-fast order, keep as-is.
    % Otherwise adapt the vectorization below to your actual pair ordering.

    % -- Geometry (vectorized) --
    % Distances TX->voxel (nx x ny x nTx)
    TXx = reshape(BP_TX_x(i,:), 1,1,[]);
    TXy = reshape(BP_TX_y(i,:), 1,1,[]);
    Rt_tx = sqrt( (Xg - TXx).^2 + (Yg - TXy).^2 + (z_offset).^2 );

    % Distances voxel->RX (nx x ny x nRx)
    RXx = reshape(BP_RX_x(i,:), 1,1,[]);
    RXy = reshape(BP_RX_y(i,:), 1,1,[]);
    Rr_rx = sqrt( (Xg - RXx).^2 + (Yg - RXy).^2 + (z_offset).^2 );

    % Combine to (nx x ny x nTx x nRx): distance sums for each (tx,rx)
    % bsxfun semantics via implicit expansion in recent MATLAB
    Rsum_txrx = Rt_tx(:,:,:,1) + Rr_rx(:,:,1,:);    % (nx,ny,nTx,nRx)

    % Expand over frequency to (nx x ny x nTx x nRx x nF)
    % Phase term k = 2*pi*f/c
    kf = reshape(2*pi*single(freq)/c, 1,1,1,1,[]);
    Rsum_full = Rsum_txrx(:,:,:,:,1) .* ones(1,1,1,1,nF,'like',Rsum_txrx); % cheap expand
    phase = exp(-1i .* kf .* Rsum_full);                                    % (nx,ny,nTx,nRx,nF)

    % Reshape to (Nobs x V) so G'*recording(:) -> Vx1
    % Order: [tx, rx, f] must match how "recording(:)" is stacked.
    % If your recording is (nTx*nRx, nF) with TX-fast, RX-slow, and then f,
    % stack as tx-fast, rx-next, f-last:
    phase = permute(phase, [3 4 5 1 2]);           % (nTx,nRx,nF,nx,ny)
    G = reshape(phase, [], V);             % (Nobs, V)

    % Accumulate image
    Im = Im + G' .* recording(:);                   % (V,1)
end

% Reshape to 2D image
Im = reshape(Im, [nx, ny]);    % xgrid rows, ygrid cols

% Magnitude and normalization
Im2d = abs(Im);
Im2d = Im2d / max(Im2d(:) + eps);

% Display magnitude (dB)
figure(20); imagesc(xgrid, ygrid, 20*log10(Im2d.'));
axis xy; xlabel('x (m)'); ylabel('y (m)');
caxis([-23 0]); colormap hot; colorbar

% Phase (radians)
figure(21); imagesc(xgrid, ygrid, angle(Im).');
axis xy; xlabel('x (m)'); ylabel('y (m)'); colorbar

% PDP from one channel (fix indexing and df)
% recording is (nTx*nRx, nF) — pick a valid row index, e.g., 250 if it exists
row_idx = min(250, size(recs,1));
df = mean(diff(freq));                     % Hz
PDP = ifft(recs(row_idx,:,1), 512, 2);     % take a scan (here i=1), FFT length 512
t = (0:numel(PDP)-1)./(df*numel(PDP));     % seconds
R = c*t/2;                                 % meters
figure; plot(R, abs(PDP)); xlabel('Range (m)'); ylabel('|PDP|');

% Optional smoothing
Im3d = imgaussfilt(Im2d, 0.7);
Im3d = Im3d / max(Im3d(:) + eps);

figure(22); imagesc(xgrid, ygrid, 20*log10(Im3d.'));
axis xy; xlabel('x (m)'); ylabel('y (m)');
caxis([-32 0]); colormap hot; colorbar
