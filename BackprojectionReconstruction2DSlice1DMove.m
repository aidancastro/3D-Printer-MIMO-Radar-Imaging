% y axis - horizontal printer move
% x axis - vetical printer move
% z axis - towards target

% make just plotting file too 
%close all; clc; clear;
disp('Reconstruction Starting'); tic;

%% ---------------- USER PARAMS ----------------------------------------
dataFile = 'XY_1D_20250724_173414.mat';
load(dataFile)   % contains variable `recs`
[~, baseName] = fileparts(dataFile);

vtrigU_ants_location;
[Xgrid,Ygrid,Zgrid]=meshgrid(xgrid,ygrid,zgrid);

src = reshape(cat(4,Xgrid,Ygrid,Zgrid),[],3);
src2 = permute(src,[3,2,4,1]);
c = physconst('lightspeed'); %(m/s)
N_freq = length(freq);
Nfft = 2^(ceil(log2(size(freq,2)))+1);

%% ---------------------- dynamic grid collapse -------------------------
gv   = {xgrid, ygrid, zgrid};           % 1: x  2: y  3: z
big  = cellfun(@numel,gv) > 10;         % “dense” axes logical mask

%‑‑ keep dense axes, collapse the sparse one to its first value
gv(~big) = cellfun(@(v){v(1)}, gv(~big), 'uni', false);
[xgrid,ygrid,zgrid] = deal(gv{:});      % overwrite originals

%‑‑ identify which dim was collapsed (needed for slicing later)
sliceDim = find(~big);                  % scalar 1|2|3

%% ---------------------- meshgrid & slice size -------------------------
[Xgrid,Ygrid,Zgrid] = meshgrid(xgrid,ygrid,zgrid);

% Map dense axes → row / col order: always rows = 1st kept axis,
%                                    cols = 2nd kept axis
denseDims      = find(big);            % two numbers, e.g. [1 2] or [1 3] …
gridCell       = {xgrid, ygrid, zgrid};
gridA          = gridCell{denseDims(1)};   % rows   (Ny)
gridB          = gridCell{denseDims(2)};   % cols   (Nx)
Ny = numel(gridA); Nx = numel(gridB);

%% ---------------------- pre‑allocate reconstruction stack -------------
nRecs   = size(recs,3);                 % how many sweeps in .mat
y_accum = zeros(Ny, Nx, nRecs, 'single');   % or 'double'
y_cart_sum = zeros(Ny, Nx, 'single');
toDB = @(M) 20*log10(abs(M)+eps) - max(20*log10(abs(M(:))+eps)); %db helper function


%% Back Projection Loop 
for i=1:nRecs
X = recs(:,:,i);
Rvec = src2-(VtrigU_ants_location + [0 0.02*(i-1) 0]);
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

%Remove resonant frequencies > =(TxRxPairs,
X = X .* (1-f_res);  
%convert to complex time domain signal
x = ifft(X,Nfft,2);
       
y_cart = reshape(H2*reshape(X,[],1),size(Xgrid));

temp_y = getSlice(y_cart, sliceDim);  % Ny × Nx matrix
y_accum(:,:,i) = temp_y;              % store this record
y_cart_sum     = y_cart_sum + temp_y; % accumulate for reconstruction

a = sprintf('Image %d/%d Proccessed',i, numel(recs(1,1,:)) ); disp(a);

end %% End of Backprojection
%% ---------- decide which plane we have --------------------------------
axisLabels = 'XYZ';            % 1:X  2:Y  3:Z

if     numel(xgrid) <= 10
        sliceName = 'YZ';  rows = ygrid; cols = zgrid;
elseif numel(ygrid) <= 10
        sliceName = 'XZ';  rows = xgrid; cols = zgrid;
else
        sliceName = 'XY';  rows = xgrid; cols = ygrid;
end

%% ---------- Improved -------------------------------------------
C_imp = y_cart_sum;     % or rssq(y_accum,3) for RSS average
plotSlice(C_imp, rows, cols, sliceName, 'Improved', baseName)

%% ---------- Unimproved (first record) ---------------------------------
C_raw = y_accum(:,:,1);
plotSlice(C_raw, rows, cols, sliceName, 'Unimproved', baseName)

%% Helper Functions
function C = getSlice(Y, sliceDim)
    switch sliceDim
        case 1,  C = squeeze(Y(1,:,:));   % collapsed Y → keep X‑Z
        case 2,  C = squeeze(Y(:,1,:));   % collapsed X → keep Y‑Z
        case 3,  C = squeeze(Y(:,:,1));   % collapsed Z → keep Y‑X
    end
end

function plotSlice(C, rows, cols, sliceName, tag, baseName)
    % Convert to dB
    Cdb = 20*log10(abs(C) + eps);

    figure
    imagesc(cols, rows, Cdb.');        % vectors → transpose
    axis xy equal tight
    colormap turbo; colorbar

    % Dynamic 5 dB limits
    lim = caxis;
    lim(1) = floor(lim(1)/5)*5;
    lim(2) = ceil(lim(2)/5)*5;
    caxis(lim)

    xlabel(sprintf('%c [m]', sliceName(2)))
    ylabel(sprintf('%c [m]', sliceName(1)))
    title(sprintf('%s Power Slice (%s)', sliceName, tag))

    outFile = sprintf('%s_%s_%s.png', baseName, sliceName, tag);
    exportgraphics(gcf, outFile, 'Resolution', 300)
end
