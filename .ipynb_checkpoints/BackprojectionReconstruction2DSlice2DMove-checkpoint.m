% y axis - horizontal printer move
% x axis - vetical printer move
% z axis - towards target

% Backprojection based off of grid center (assume NxN square grid)
%close all; clc; clear;

%% ---------------- USER PARAMS ----------------------------------------

close all; clear; disp('Reconstruction Starting'); tic;

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'YZ_2D_20250731_161912.mat';
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
Ny = numel(xgrid);
Nx = numel(ygrid);
Nz = numel(zgrid);

%% pre‑allocate reconstruction stack -------------
nRecs   = size(y_cart,3);                 % how many sweeps in .mat
y_accum = zeros(Ny, Nx, Nz, 'single');   % or 'double'
y_cart_sum = zeros(Ny, Nx, 'single');
toDB = @(M) 20*log10(abs(M)+eps) - max(20*log10(abs(M(:))+eps)); %db helper function


%% Back Projection Loop --------------------------
for i=1:nRecs
X = y_cart(:,:,i);
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
padsig = padsig(:).';                 % ensure row
k = (lnconv-1)/2;
padsig = padarray(padsig, [0 k], 'symmetric', 'both');
padsig = conv(padsig, c2, 'valid');

f_res = padsig > thresh;

%Remove resonant frequencies
X = X .* (1-f_res);  
%convert to complex time domain signal
x = ifft(X,Nfft,2);
       
y_cart = reshape(H2.*reshape(X,[],1),size(Xgrid));

% pull out the slice (will be linear if one dim was collapsed)
temp_y = squeeze(y_cart);    % → [Ny×Nx] matrix every time

% now safe to accumulate


y_cart_sum = y_cart_sum + temp_y; % accumulate for reconstruction

[m,~] = cfar2d_localmax(temp_y, 1, 6, 1e-4, 1);
y_accum(:,:,i) = temp_y .* m; 
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
%plotSlice(Norm_C_imp, rows, cols, sliceName, 'Improved', baseName)

%% Unimproved (first record) ---------------------------------
C_raw = y_accum(:,:,ceil(nRecs/2) ); %return grid center image
Norm_C_raw = C_raw ./ max(abs(C_raw(:)));
%plotSlice(Norm_C_raw, rows, cols, sliceName, 'Unimproved', baseName)

% rows, cols, sliceName already set by your if/elseif block
for i = 1:nRecs:1
   Temp_C = y_accum(:,:,i);
   Normalize_C = Temp_C ./ max(abs(Temp_C(:)));
   C_dB = 20*log10(abs(Normalize_C) + eps); 
   C_all(:,:,i) = Normalize_C;
end
hFig = plotSlicesGridByPlane(stack_db, rows, cols, sliceName, ...
        struct('Title',"All slices", 'Scale',"db"));


%% Helper Functions

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
function [mask,T] = cfar2d_localmax(P_lin, g, t, pfa, r)
% P_lin: 2-D power (linear). g: guard half-width. t: training half-width.
% pfa: desired false-alarm rate (e.g., 1e-4). r: local-max radius.
P = double(P_lin);

Kall   = ones(2*(g+t)+1);
Kguard = ones(2*g+1);

sumAll   = conv2(P, Kall,   'same');
sumGuard = conv2(P, Kguard, 'same');
Ntrain   = numel(Kall) - numel(Kguard);
noiseHat = max(sumAll - sumGuard, 0) / max(Ntrain,1);

alpha = Ntrain * (pfa^(-1/Ntrain) - 1);   % CA-CFAR scale
T     = alpha .* noiseHat;

isLM  = imdilate(P, true(2*r+1)) == P;    % non-maximum suppression
mask  = (P > T) & isLM;
end
function hFig = plotSlicesGridByPlane(stack, rows, cols, sliceName, opts)
% stack: [numel(rows) x numel(cols) x N]  (each slice already oriented to rows x cols)
% rows, cols: coordinate vectors used in your switch above
% sliceName: 'XY','YX','XZ','ZX','YZ','ZY'  (first = vertical, second = horizontal)
% opts (all optional): opts.Scale="db"|"linear", opts.CLim=[lo hi], opts.Title, ...
%                      opts.Order (permute slices), opts.Colormap="parula", opts.SaveAs

arguments
    stack {mustBeNumeric, mustBeNonempty}
    rows (:,1) double
    cols (:,1) double
    sliceName (1,2) char
    opts.Scale (1,1) string = "db"
    opts.CLim (1,2) double = []
    opts.Title (1,1) string = ""
    opts.Order double = []
    opts.Colormap (1,1) string = "parula"
    opts.SaveAs (1,1) string = ""
end

[ny,nx,N] = size(stack);
if numel(rows)~=ny || numel(cols)~=nx
    error('Size mismatch: stack is [%d x %d x %d], rows=%d, cols=%d.', ny,nx,N,numel(rows),numel(cols));
end

% Convert to dB for display if needed
if opts.Scale=="linear"
    stack_db = 10*log10(stack + eps);
else
    stack_db = stack;
end

% Global color limits (robust to outliers)
if isempty(opts.CLim)
    v = stack_db(isfinite(stack_db));
    clim = [prctile(v,2) prctile(v,98)];
else
    clim = opts.CLim;
end

% Slice order
if isempty(opts.Order), order = 1:N; else, order = opts.Order(:)'; end

% Grid size
ncols = ceil(sqrt(N));
nrows = ceil(N/ncols);

% Labels from sliceName: first char = vertical (rows), second = horizontal (cols)
ylab = sprintf('%s [m]', sliceName(1));
xlab = sprintf('%s [m]', sliceName(2));

% Figure + layout
hFig = figure('Color','w','Position',[100 100 1200 900]);
tl = tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','compact');
if strlength(opts.Title)>0, title(tl, opts.Title); end
colormap(opts.Colormap);

for k = 1:N
    kk = order(k);
    ax = nexttile;
    imagesc(cols, rows, stack_db(:,:,kk));     % <-- rows on vertical, cols on horizontal
    axis image; set(ax,'YDir','normal'); caxis(clim);
    title(sprintf('Slice %d',kk),'FontSize',8);

    % Clean inner ticks
    if mod(k-1,ncols)~=0, ax.YTickLabel = []; end
    if k <= (nrows-1)*ncols, ax.XTickLabel = []; end

    % Edge labels only
    if mod(k-1,ncols)==0, ylabel(ylab); end
    if k > (nrows-1)*ncols, xlabel(xlab); end
end

cb = colorbar; cb.Layout.Tile = 'east'; ylabel(cb,'Power [dB]');

if strlength(opts.SaveAs)>0
    exportgraphics(hFig, opts.SaveAs, 'Resolution', 300);
end
end

toc
