function Backprojection2D()
% Backprojection2D — 2D backprojection from an NxN scan after collapsing
% the sparsest axis. Allocation is adaptive to the remaining two axes.
%
% Requires in the loaded .mat or on path:
%   recs (freq x TxRx x nRecs), freq, xgrid, ygrid, zgrid,
%   TxRxPairs, VtrigU_ants_location (function/vars), printer_offsets,
%   RadiationPattern(theta,phi)
%
% Outputs (left in workspace):
%   y_accum [numel(gridA) x numel(gridB) x nRecs]
%   y_cart_sum [numel(gridA) x numel(gridB)]
%   rows, cols, sliceName, axisLabels

disp('Reconstruction Starting'); tic;

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'YZ_2D_20250811_204940.mat';
load(dataFile);                          % expects recs, freq, etc.
[~, baseName] = fileparts(dataFile);     %#ok<NASGU>

% Antenna geometry / grids / pairs / patterns
vtrigU_ants_location;    % should define xgrid,ygrid,zgrid, TxRxPairs, VtrigU_ants_location

c      = physconst('lightspeed');
N_freq = numel(freq(:));
lambda = c ./ freq(:);                   % [N_freq x 1]
Nfft   = 2^(ceil(log2(N_freq))+1);       %#ok<NASGU> % kept if you later need IFFT

%% -------- dynamic grid collapse (collapse the smallest axis) ----------
lens = [numel(xgrid), numel(ygrid), numel(zgrid)];
[~, sliceDim] = min(lens);               % 1=X, 2=Y, 3=Z

gv = {xgrid, ygrid, zgrid};
gv{sliceDim} = gv{sliceDim}(1);          % collapse sparsest axis to first value
[xgrid, ygrid, zgrid] = deal(gv{:});

% Rebuild grids after collapse
[Xgrid, Ygrid, Zgrid] = meshgrid(xgrid, ygrid, zgrid);
V = numel(Xgrid);
vox = [Xgrid(:), Ygrid(:), Zgrid(:)];    % [V,3] voxel list

% The two dense axes define rows/cols
denseDims = setdiff(1:3, sliceDim);
gridA     = gv{denseDims(1)};            % rows
gridB     = gv{denseDims(2)};            % cols

%% ---------------- adaptive allocation ----------------
nRecs      = size(recs, 3);
y_accum    = zeros(numel(gridA), numel(gridB), nRecs, 'single');
y_cart_sum = zeros(numel(gridA), numel(gridB), 'single');

%% ------------------- Back Projection Loop ----------------------------
% constants for propagation
RCS  = 1;                                % m^2
csf  = sqrt(RCS) .* lambda ./ ((4*pi).^(3/2));   % [N_freq x 1]

for i = 1:nRecs
    X = recs(:, :, i);                   % [N_freq x (TxRx)]
    % Antenna locations for this printer pose
    ant_xyz = VtrigU_ants_location + [printer_offsets(i,1), printer_offsets(i,2), 0];

    %% -- Identify resonant frequencies (your method) --
    thresh = 3;
    df = abs(freq(2) - freq(1));         % Hz (assumes uniform)
    lnconv = min( ...
        max(floor(N_freq/8)*2+1, max(3, floor(50/df)*2+1)), ...
        floor(3*N_freq/8)*2+1 );
    c2 = -ones(lnconv,1)/(lnconv-1);
    c2((lnconv+1)/2) = 1;

    padsig = 20*log10(rssq(X,1));
    padsig = [padsig((lnconv-1)/2:-1:1), padsig, padsig(end:-1:end-(lnconv-1)/2+1)];
    padsig = conv(padsig, c2, 'valid');
    f_res  = padsig > thresh;            % logical mask (1 x TxRx)
    X      = X .* (1 - f_res);           % zero-out resonant freqs per channel

    %% -- Build steering/propagation matrix H2: [V  x (TxRx*N_freq)] --
    H2 = zeros(V, numel(TxRxPairs)*N_freq, 'like', 1+1j);
    col = 1;

    for ii = 1:size(TxRxPairs,1)
        tx = TxRxPairs(ii,1);
        rx = TxRxPairs(ii,2);

        % Distances and angles from tx/rx to every voxel
        vtx = vox - ant_xyz(tx,:);  Rtx = vecnorm(vtx,2,2);
        vrx = vox - ant_xyz(rx,:);  Rrx = vecnorm(vrx,2,2);

        Rtheta_tx = atan2(vecnorm(vtx(:,1:2),2,2), vtx(:,3));
        Rphi_tx   = atan2(vtx(:,2), vtx(:,1));
        Rtheta_rx = atan2(vecnorm(vrx(:,1:2),2,2), vrx(:,3));
        Rphi_rx   = atan2(vrx(:,2), vrx(:,1));

        % Pattern magnitudes (user function)
        Smag_tx = RadiationPattern(Rtheta_tx, Rphi_tx) ./ max(Rtx, eps);
        Smag_rx = RadiationPattern(Rtheta_rx, Rphi_rx) ./ max(Rrx, eps);

        % Two-way phase for each frequency: 2*pi*(Rtx+Rrx)/lambda
        % Expand to [V x N_freq]
        Sphase = 2*pi * (Rtx + Rrx) .* (1 ./ lambda.');      % broadcast over freq

        % Combined scalar factor per freq
        gain = 1 ./ (csf.' .* (Smag_tx .* Smag_rx));         % [1 x N_freq]

        % Complex propagation for this Tx/Rx pair
        Hpair = gain .* exp(-1j * Sphase);                   % [V x N_freq]

        % Place into big H2
        H2(:, col:col+N_freq-1) = Hpair;
        col = col + N_freq;
    end

    %% -- Backproject --
    % Flatten data to match H2 columns: [(TxRx*N_freq) x 1]
    Xcol = reshape(X.', [], 1);                 % pair-major then freq

    y_vol = reshape(H2 * Xcol, size(Xgrid));    % [Ny x Nx x Nz] with one dim == 1

    %% -- Extract 2D slice with correct orientation for rows/cols --
    switch sliceDim
        case 1      % collapsed X  -> YZ plane (rows=y, cols=z)
            temp_y = squeeze(y_vol(:,1,:));                 % [numel(ygrid) x numel(zgrid)]
        case 2      % collapsed Y  -> XZ plane (rows=x, cols=z)
            temp_y = squeeze(y_vol(1,:,:));                 % [numel(xgrid) x numel(zgrid)]
        case 3      % collapsed Z  -> XY plane (rows=x, cols=y)
            temp_y = squeeze(y_vol(:,:,1)).';               % transpose -> rows=x, cols=y
    end

    % Accumulate
    y_accum(:,:,i) = single(temp_y);
    y_cart_sum     = y_cart_sum + single(temp_y);

    fprintf('Image %d/%d Processed\n', i, nRecs);
end

%% ---------------- plane labels for downstream plotting ---------------
axisLabels = 'XYZ';   % 1:X  2:Y  3:Z
switch sliceDim
    case 1, sliceName = 'YZ'; rows = ygrid; cols = zgrid;
    case 2, sliceName = 'XZ'; rows = xgrid; cols = zgrid;
    case 3, sliceName = 'XY'; rows = xgrid; cols = ygrid;
end

toc
end
