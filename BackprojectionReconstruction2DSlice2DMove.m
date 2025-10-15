% y axis - horizontal printer move
% x axis - vertical printer move
% z axis - towards target

clear;
close all; disp('Reconstruction Starting'); tic;


%load('Calibration_XY_2D_20251005_170345.mat');
%calibration_recs = recs(:,:,:);

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'YZ_2D_20251014_215342.mat';
load(dataFile)   % expects: recs, freq, xgrid, ygrid, zgrid, TxRxPairs, VtrigU_ants_location, RadiationPattern
[~, baseName] = fileparts(dataFile);

%% ---- Antenna locations (x=vertical, y=horizontal, z=0) ----
% File created from Python: contains 20×1 vectors TX_x, TX_y, RX_x, RX_y
S = load('AntennaLocations.mat');   % adjust path if needed
TX_x = double(S.TX_x(:));           % TX 1..20 (top→bottom)
TX_y = double(S.TX_y(:));           % constant y for TX
RX_x = double(S.RX_x(:));           % constant x for RX
RX_y = double(S.RX_y(:));           % RX 21..40 (left→right)

% Hardware order on VK-74 is RX(1..20), then TX(21..40)
ant_x = [TX_x; RX_x];               % (40×1)
ant_y = [TX_y; RX_y];               % (40×1)
ant_z = zeros(40,1);                % antenna plane (z=0)

% beam-steering math expects [y x z]; keep your original packing
Vant = [ant_y, ant_x, ant_z];        % (40×3)  [y x z]

%%% <<< FIX #1: ensure antenna units are meters (file is often in mm)
if max(abs(Vant(:))) > 1            % heuristic: values look like mm if > 1
    Vant = Vant * 1e-3;             % mm -> m
end

% -[Xgrid, Ygrid, Zgrid] = meshgrid(xgrid, ygrid, zgrid);
% -src  = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);  % (V,3) [X Y Z]
% -src2 = permute(src,[3,2,4,1]);

[Xgrid, Ygrid, Zgrid] = meshgrid(xgrid, ygrid, zgrid);
src  = [Xgrid(:), Ygrid(:), Zgrid(:)];                 % (V,3) [X Y Z]
src2 = reshape(src.', [1,3,1,numel(src)/3]);           % (1×3×1×V)


c     = physconst('lightspeed'); %(m/s)
N_freq = numel(freq);
Nfft   = 2^(ceil(log2(N_freq))+1);

%% -------- collapse the thinnest axis, then (re)build grids & src --------
lens = [numel(xgrid), numel(ygrid), numel(zgrid)];
[~, sliceDim] = min(lens);

gv = {xgrid, ygrid, zgrid};
gv{sliceDim} = gv{sliceDim}(1);           % collapse to length 1
[xgrid, ygrid, zgrid] = deal(gv{:});

% Decide plane + rows/cols consistent with [Y X Z] packing
if     sliceDim == 3          % Z collapsed -> XY slice
    sliceName = 'XY'; rows = ygrid; cols = xgrid;
elseif sliceDim == 2          % Y collapsed -> XZ slice
    sliceName = 'XZ'; rows = xgrid; cols = zgrid;
else                          % X collapsed -> YZ slice
    sliceName = 'YZ'; rows = ygrid; cols = zgrid;
end
Ny = numel(rows); Nx = numel(cols);

% Rebuild grids and voxel list AFTER collapse
[Ygrid, Xgrid, Zgrid] = meshgrid(ygrid, xgrid, zgrid);   % note order
src  = [Ygrid(:), Xgrid(:), Zgrid(:)];                   % (V,3) [Y X Z]
src2 = reshape(src.', [1,3,1,numel(src)/3]);             % (1×3×1×V)

% Pre-allocate using the 2-D slice size
nRecs      = size(recs,3);
y_accum    = zeros(Ny, Nx, nRecs, 'single');
y_cart_sum = zeros(Ny, Nx, 'single');

toDB = @(M) 20*log10(abs(M)+eps) - max(20*log10(abs(M(:))+eps));


% 5x5 snake printer offsets ([y x]) in METERS

% printer_offset_y = [0 step 2*step 3*step 4*step 4*step -3*step -2*step -step 0 0 -step -2*step -3*step -4*step -4*step -3*step -2*step -step 0 0 -step -2*step -3*step -4*step];
% printer_offset_x = [0 0 0 0 0 -step -step -step -step -step -2*step -2*step -2*step -2*step -2*step -3*step -3*step -3*step -3*step -3*step -4*step -4*step -4*step -4*step -4*step];
% printer_offsets_yx  = [printer_offset_y(:), printer_offset_x(:)];
% --- Snake pattern offsets for 5×5 scan ---
N = 5;

step = 0.02;  % 20 mm = 0.02 m

printer_offset_y = zeros(N^2,1);  % stage motion along Y axis
printer_offset_x = zeros(N^2,1);  % stage motion along X axis

for k = 1:N^2
    row = floor((k-1)/N);  % index of X step (0–4)
    col = mod(k-1, N);     % index of Y step (0–4)

    if mod(row,2) == 0
        % even row: move in +Y direction
        y =  col;
    else
        % odd row: move in -Y direction
        y = (N-1) - col;
    end

    printer_offset_x(k) =  row * step;  % +X each new row
    printer_offset_y(k) =  y   * step;  % ±Y within each row
end

% combine in [Y X] order (physical stage positions)
printer_offsets_yx = [printer_offset_y, printer_offset_x];


% grid-center correction (these were in mm in your notes)
firstscan_off = [-(gridCenter(1,1)-2*xstep), gridCenter(1,2)-2*xstep, 0] * 1e-3; % m
Vant = Vant + firstscan_off;

y_ref = [];   % reference complex image for phase locking

%% Back Projection Loop
for i = 1:nRecs
    X = recs(:,:,i);  % (nPairs, nFreq)
%    X = X - calibration_recs(:,:,i);

       % pack antennas once for 4-D implicit expansion
    K = size(Vant,1);
    Vant43 = reshape(Vant.', [1, 3, K, 1]);     % (1×3×K×1)  [Y X Z]
    
    % per-scan printer offset (still in [Y X 0] to match Vant)
    vox_shift_vec = reshape([printer_offsets_yx(i,1), ...
                             printer_offsets_yx(i,2), 0], [1,3,1,1]);  % (1×3×1×1)
    
    % shift voxels opposite platform motion, then subtract antenna coords
    Rvec = (src2 - vox_shift_vec) - Vant43;     % (1×3×K×V)

    
    % Now make distance/angles with consistent dims
    Rmag  = sqrt(sum(Rvec.^2, 2));             % (1×1×K×V)
    Rmag  = permute(Rmag, [3 2 1 4]);          % -> (K×1×1×V)
    
    Rxy   = sqrt(sum(Rvec(:,1:2,:,:).^2, 2));  % (1×1×K×V)
    Rxy   = permute(Rxy,  [3 2 1 4]);          % -> (K×1×1×V)
    
    Ry    = permute(Rvec(:,2,:,:), [3 2 1 4]); % (K×1×1×V)
    Rx    = permute(Rvec(:,1,:,:), [3 2 1 4]); % (K×1×1×V)
    Rz    = permute(Rvec(:,3,:,:), [3 2 1 4]); % (K×1×1×V)
    
    Rtheta = atan2(Rxy, Rz);                   % (K×1×1×V)
    Rphi   = atan2(Ry, Rx);                    % (K×1×1×V)

    % --- steering terms (keep your original forms) ---
    nPairs = size(TxRxPairs,1);
    nF     = numel(freq);
    V      = size(src2,4);            % voxels

    Sphase = 2*pi*Rmag .* reshape(freq(:).', [1 nF 1 1]) / c;  % (K×nF×1×V)
    Sphase = permute(Sphase, [1 2 3 4]);                       % unchanged

    Smag = 10^(5.8/20) * RadiationPattern(Rtheta, Rphi) ./ max(Rmag, eps('like',Rmag)); % (K×1×1×V)

    lambda = c ./ freq(:).';
    csf    = reshape((sqrt(1) .* lambda) ./ ((4*pi).^(3/2)), [1 nF 1 1]);

    H2 = complex(zeros(nPairs, nF, 1, V, 'like', X));          % (nPairs×nF×1×V)
    for ii = 1:nPairs
        tx = TxRxPairs(ii,1);  rx = TxRxPairs(ii,2);
        H2(ii,:,:,:) = 1 ./ ( csf .* Smag(tx,1,1,:) .* Smag(rx,1,1,:) .* exp(-1j*( Sphase(tx,:,:,:) + Sphase(rx,:,:,:) )) );
    end
    H2 = reshape(permute(H2, [4 1 2 3]), V, nPairs*nF);        % (V × nPairs*nF)


    %resonance removal
    thresh = 3; Nf = size(X,2);
    if Nf >= 3
        df = (Nf>=2) * abs(median(diff(freq))) + (Nf<2)*1;
        cand = [2*floor(Nf/8)+1, 2*floor(50/df)+1, 2*floor(3*Nf/8)+1];
        lnconv = min([cand, Nf - 1 - mod(Nf-1,2)]);
        lnconv = max(3, lnconv); if mod(lnconv,2)==0, lnconv = lnconv - 1; end
        Lh = (lnconv-1)/2;
        c2 = -ones(lnconv,1,'like',real(X))/(lnconv-1); c2(Lh+1) = 1;
        mags  = 20*log10(max(rssq(X,1), eps('like',X)));
        mags  = double([mags(Lh:-1:1), mags, mags(end:-1:end-Lh+1)]);
        trace = conv(mags, c2, 'valid');      % length == Nf
        f_res = trace > thresh;                % 1×Nf logical
        X(:, f_res) = 0;
    end

    % ---- accumulate ----
    y_cart = reshape(H2 * X(:), [numel(rows), numel(cols)]);

    % ---- NEW: global phase alignment (one scalar per scan) ----
    if isempty(y_ref)
        % First scan becomes reference; or use the running average
        y_ref = y_cart;
    else
        % Estimate constant phase offset that best aligns scan i to the reference
        phi = angle( sum( conj(y_ref(:)) .* y_cart(:) ) );   % [-pi, pi]
        y_cart = y_cart * exp(-1j*phi);
    end
    % -----------------------------------------------------------

y_accum(:,:,i) = squeeze(y_cart);
y_cart_sum     = y_cart_sum + y_accum(:,:,i);


    fprintf('Image %d/%d processed\n', i, nRecs);
end

%% --------- normalization + sums (unchanged) ----------
[Ssum, Snc, y_norm] = norm_sum_simple(y_accum);

C_imp      = Ssum;
Norm_C_imp = C_imp ./ max(abs(C_imp(:)) + eps);
plotSlice(Norm_C_imp, rows, cols, sliceName, 'Improved (coherent, normalized)', baseName)

midIdx     = 1;
C_raw      = y_norm(:,:,midIdx);
Norm_C_raw = C_raw ./ max(abs(C_raw(:)) + eps);
plotSlice(Norm_C_raw, rows, cols, sliceName, 'Unimproved (single scan, normalized)', baseName)

toc
outFile = sprintf('reconstructed_%s_%s.mat', sliceName, baseName);
%save(outFile, 'y_accum', 'y_norm', 'Ssum', 'Snc', 'xgrid','ygrid','zgrid','Xgrid','Ygrid','Zgrid');

% Use complex images before display scaling
Img_imp = C_imp;               % improved (coherent sum)
Img_raw = C_raw;               % original (single scan), from y_norm(:,:,midIdx)

% Work with magnitudes (power later)
Mag_imp = abs(Img_imp);
Mag_raw = abs(Img_raw);

% --- Find peak on the improved image and define ROIs ---
% Peak location (row, col)
[~, linMax] = max(Mag_imp(:));
[peak_r, peak_c] = ind2sub(size(Mag_imp), linMax);

% Set ROI radii in pixels (tune r_sig/r_in/r_out to your resolution)
r_sig = 2;         % signal disk radius (pixels)
r_in  = 5;         % inner radius of noise annulus
r_out = 12;        % outer radius of noise annulus

% Build distance map in pixel units
[CC, RR] = meshgrid(1:size(Mag_imp,2), 1:size(Mag_imp,1));
D = hypot(RR - peak_r, CC - peak_c);

sigMask  = (D <= r_sig);
noiseMask = (D >= r_in) & (D <= r_out);

% Safety: ensure we have noise pixels
if nnz(noiseMask) < 50
    warning('Noise annulus too small—expand r_out or shrink r_in.');
end

% --- Helper to compute SNR (dB) for any image given fixed ROIs ---
snr_db = @(Img) 10*log10( mean( abs(Img(sigMask)).^2 ) / ...
                          mean( abs(Img(noiseMask)).^2 + eps ) );

SNR_imp_dB = snr_db(Img_imp);
SNR_raw_dB = snr_db(Img_raw);

fprintf('SNR (improved, coherent sum): %.2f dB\n', SNR_imp_dB);
fprintf('SNR (original, single scan) : %.2f dB\n', SNR_raw_dB);

%% ----------------- Helpers (unchanged) -----------------
function [Ssum, Snc, Snorm] = norm_sum_simple(S)
[Ny,Nx,N] = size(S);
Snorm = zeros(Ny,Nx,N,'like',S);
sc = zeros(N,1,'like',real(S));
epsm = 1e-12;
for i = 1:N
    A = S(:,:,i);
    mag = abs(A);
    cap = prctile(mag(:), 99);
    A   = A .* min(1, cap./max(mag,epsm));
    p90  = prctile(abs(A(:)), 90);
    pool = abs(A(abs(A) <= p90)); if isempty(pool), pool = abs(A(:)); end
    sc(i) = median(pool) + epsm;
    Snorm(:,:,i) = A;
end
ref = median(sc); g = ref ./ sc;
for i = 1:N, Snorm(:,:,i) = Snorm(:,:,i) .* g(i); end
Ssum = sum(Snorm,3);
Snc  = sqrt(sum(abs(Snorm).^2,3));
end

function plotSlice(C, rows, cols, sliceName, tag, ~)

    figure;
    C2  = squeeze(C);  assert(ndims(C2)==2, 'plotSlice: expected 2-D input');
    Cdb = 20*log10(abs(C2) + eps); 
    Cdb = max(Cdb - max(Cdb(:)), -40);

    switch sliceName
        case 'XY'
            x = cols; y = rows; xlab = 'X [m]'; ylab = 'Y [m]';
        case 'YZ'
            x = cols; y = rows; xlab = 'Z [m]'; ylab = 'Y [m]';
        case 'XZ'
            x = cols; y = rows; xlab = 'Z [m]'; ylab = 'X [m]';
        otherwise
            x = cols; y = rows; xlab = 'cols';   ylab = 'rows';
    end

    imagesc(y, x, Cdb);          % <-- no permute; x first, then y
    axis xy equal tight; 
    colormap turbo; colorbar; clim([-3 0])
    xlabel(ylab); ylabel(xlab);
    title(sprintf('%s Power Slice — %s', sliceName, tag))
end

