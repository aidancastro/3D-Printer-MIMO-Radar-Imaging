%% --- Antenna‑location visualiser --------------------------------------
% VtrigU_ants_location must exist in the workspace
vtrigU_ants_location;
A = VtrigU_ants_location;          % (m)  N×3  [x y z]

txIdx = 1:20;                      % first 20 are Tx, rest Rx  (adapt if needed)
rxIdx = 21:40;

figure('Name','Antenna positions, VtrigU board'); clf;

% --- 3‑D scatter --------------------------------------------------------
subplot(1,2,1); hold on; grid on; axis equal;
scatter3(A(txIdx,1),A(txIdx,2),A(txIdx,3),60,'r','filled');
scatter3(A(rxIdx,1),A(rxIdx,2),A(rxIdx,3),60,'b','filled');

% coordinate‑frame axes for easy orientation check
L = 0.04;                                            % 4 cm arrow length
quiver3(0,0,0, L,0,0,'k','LineWidth',1.5,'MaxHeadSize',0.5);  % +x
quiver3(0,0,0, 0,L,0,'k','LineWidth',1.5,'MaxHeadSize',0.5);  % +y
quiver3(0,0,0, 0,0,L,'k','LineWidth',1.5,'MaxHeadSize',0.5);  % +z
text(L,0,0,' +x'); text(0,L,0,' +y'); text(0,0,L,' +z');

xlabel('x  (m)'); ylabel('y  (m)'); zlabel('z  (m)');
title('3‑D view'); view(35,25);
legend({'Tx 1–20','Rx 21–40'},'Location','northeast');

% --- top view (x–y plane) ----------------------------------------------
subplot(1,2,2); hold on; grid on; axis equal;
plot(A(txIdx,1),A(txIdx,2),'r.','MarkerSize',20);
plot(A(rxIdx,1),A(rxIdx,2),'b.','MarkerSize',20);
xlabel('x  (m)'); ylabel('y  (m)');
title('Top‑down (z eliminated)'); view(2);   % 2‑D view

% -----------------------------------------------------------------------
