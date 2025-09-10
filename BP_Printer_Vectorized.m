clear all; close all;

%% Load
load("Gun_Scan_20250908_201059.mat");      % recs (nPairs,nF,nRecs), freq, xgrid, ygrid, printer_offsets(:,[yoff xoff])
load("AntennaLocations.mat");              % TX_x, TX_y, RX_x, RX_y

c = 299792458;
z_offset = single(0.50);

%% Sizes
nRecs = size(recs,3);
nF    = numel(freq);
nTx   = numel(TX_x);
nRx   = numel(RX_y);
nx    = numel(xgrid);
ny    = numel(ygrid);
V     = nx*ny;

%% Voxel grid
[Xg, Yg] = ndgrid(single(xgrid), single(ygrid));
xg = Xg(:).';      % (1,V)
yg = Yg(:).';      % (1,V)

%% Pair mapping: p = (rx-1)*nTx + tx  (TX fastest)
[TXi, RXi] = ndgrid(1:nTx, 1:nRx);
tx_of_p = TXi(:);               % (nPairs,1)
rx_of_p = RXi(:);               % (nPairs,1)
nPairs  = numel(tx_of_p);

%% Frequency factor
freq = single(freq(:));         % (nF,1)
wfac = (-1j * 2*pi / c) * freq; % (nF,1)

%% --- Pre-shift antenna positions for all scans (do this BEFORE the loop) ---
% printer_offsets = [yoff, xoff] for each scan
yoff = single(printer_offsets(:,1));  % (nRecs,1)
xoff = single(printer_offsets(:,2));  % (nRecs,1)

TXx0 = single(TX_x(:)).';   % (1,nTx)
TXy0 = single(TX_y(:)).';   % (1,nTx)
RXx0 = single(RX_x(:)).';   % (1,nRx)
RXy0 = single(RX_y(:)).';   % (1,nRx)

BP_TX_x = repmat(TXx0, nRecs, 1) - xoff;   % (nRecs,nTx)
BP_TX_y = repmat(TXy0, nRecs, 1) - yoff;   % (nRecs,nTx)
BP_RX_x = repmat(RXx0, nRecs, 1) - xoff;   % (nRecs,nRx)
BP_RX_y = repmat(RXy0, nRecs, 1) - yoff;   % (nRecs,nRx)
% ---------------------------------------------------------------------------

%% Accumulators
Im      = complex(zeros(V,1,'single'));          % coherent sum (global)
Im_all  = complex(zeros(nx, ny, nRecs, 'single'));  % per-scan images

for i = 1:nRecs
    % Use pre-shifted positions for scan i (row vectors -> column)
    TXx = BP_TX_x(i,:).';   % (nTx,1)
    TXy = BP_TX_y(i,:).';   % (nTx,1)
    RXx = BP_RX_x(i,:).';   % (nRx,1)
    RXy = BP_RX_y(i,:).';   % (nRx,1)

    % Ranges to all voxels for this scan
    RT = sqrt( (TXx - xg).^2 + (TXy - yg).^2 + z_offset.^2 );  % (nTx,V)
    RR = sqrt( (RXx - xg).^2 + (RXy - yg).^2 + z_offset.^2 );  % (nRx,V)

    % Spectrum for this scan
    S = single(recs(:,:,i));   % (nPairs,nF)
    Im_i = complex(zeros(V,1,'single'));

    % Accumulate over pairs (reuse RT/RR)
    for p = 1:nPairs
        L = RT(tx_of_p(p), :) + RR(rx_of_p(p), :);   % (1,V)
        E = exp(wfac * L);                           % (nF,V)

        Im   = Im   + (E.' * S(p,:).');             % (V,1) coherent global
        Im_i = Im_i + (E.' * S(p,:).');             % (V,1) this scan
    end

    Im_all(:,:,i) = reshape(Im_i, nx, ny);
    fprintf("%d of %d scan processed\n", i, nRecs);
end

% --- Plot all 25 individual scans (x vertical, y horizontal) ---

% --- (Optional) coherent-sum image with same orientation ---
Im_sum = reshape(Im, nx, ny);
Im_sum_mag = abs(Im_sum);
Im_sum_mag = Im_sum_mag ./ (max(Im_sum_mag(:)) + eps('single'));
figure('Name','Coherent Sum','NumberTitle','off');
imagesc(ygrid, xgrid, 20*log10(max(single(Im_sum_mag), eps('single'))));
axis xy; xlabel('y (m)'); ylabel('x (m)');
caxis([-23 0]); colormap hot; colorbar;
title('Coherent Sum (dB)');

%% ---------------- helper ----------------
function plot_all_scans(Im_all, xgrid, ygrid, clim_db)
% Im_all: (nx,ny,nRecs) complex images (one per scan)
% xgrid: vertical axis (rows)
% ygrid: horizontal axis (cols)
% clim_db: [min max] dB range, e.g., [-23 0]

    if nargin < 4, clim_db = [-23 0]; end
    [nx, ny, nRecs] = size(Im_all);

    % Global magnitude normalization (consistent across all tiles)
    mag = abs(Im_all);
    mag = mag ./ (max(mag(:)) + eps('single'));      % normalized to 1
    db  = 20*log10(max(single(mag), eps('single'))); % avoid log(0), keep single

    % Layout
    if nRecs == 25
        nrows = 5; ncols = 5;
    else
        nrows = ceil(sqrt(nRecs));
        ncols = ceil(nRecs/nrows);
    end

    f = figure('Name','Individual BP Images','NumberTitle','off');
    tl = tiledlayout(f, nrows, ncols, 'Padding','compact','TileSpacing','compact');

    for i = 1:nRecs
        ax = nexttile(tl);
        imagesc(xgrid, ygrid, db(:,:,i));   % X=horiz=y, Y=vert=x
        axis(ax, 'xy', 'tight');
        xlabel('y (m)'); ylabel('x (m)');
        title(sprintf('Scan %d', i), 'FontSize', 9);
        colormap(ax, hot);
        caxis(ax, clim_db);
    end

    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Intensity (dB)';
end

plot_all_scans(Im_all, xgrid, ygrid, [-23 0]);
