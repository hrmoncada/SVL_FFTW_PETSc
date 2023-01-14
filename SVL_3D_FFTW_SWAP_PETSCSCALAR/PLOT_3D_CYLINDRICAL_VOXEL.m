%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables
clf         % clear functions

% Grid size
Sx = 1; % x-axis unit cell grid size
Sy = 1; % y-axis unit cell grid size
Sz = 1; % y-axis unit cell grid size

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Loading Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('OUTPUT_UNIT_CELL_REAL')
load ('OUTPUT_FFTW_REAL')
load ('OUTPUT_FFTW_IMAG')
load ('OUTPUT_INV_FFTW_REAL')
load ('OUTPUT_INV_FFTW_IMAG')
load ('OUTPUT_SWAP_REAL')
load ('OUTPUT_SWAP_IMAG')
load ('OUTPUT_TRUNC_FFTW_REAL')
load ('OUTPUT_TRUNC_FFTW_IMAG')
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Length Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
M1_rows = length(OUTPUT_UNIT_CELL_REAL(:,1));       % Unit Cell matrix    	(Nx x Ny x Nz)
M1_cols = length(OUTPUT_UNIT_CELL_REAL(1,:));       % Unit Cell matrix    	(Nx x Ny x Nz)

M2_rows = length(OUTPUT_FFTW_REAL(:,1));       	    % FFTW REAL matrix    	(Nx x Ny x Nz)
M2_cols = length(OUTPUT_FFTW_REAL(1,:));            % FFTW REAL matrix    	(Nx x Ny x Nz)

M3_rows = length(OUTPUT_INV_FFTW_REAL(:,1));        % FFTW REAL matrix    	(Nx x Ny x Nz)
M3_cols = length(OUTPUT_INV_FFTW_REAL(1,:));        % FFTW REAL matrix    	(Nx x Ny x Nz)

M4_rows = length(OUTPUT_SWAP_REAL(:,1));            % FFTW REAL matrix    	(Nx x Ny x Nz)
M4_cols = length(OUTPUT_SWAP_REAL(1,:));            % FFTW REAL matrix    	(Nx x Ny x Nz)

M5_rows = length(OUTPUT_TRUNC_FFTW_REAL(:,1));      % TRUC FFTW REAL matrix     (NM x NN x NP)
M5_cols = length(OUTPUT_TRUNC_FFTW_REAL(1,:));      % TRUC FFTW REAL matrix     (NM x NN x NP)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Array Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
for j = 1 : M1_cols   U1(:,:,j) = OUTPUT_UNIT_CELL_REAL((j-1)*M1_cols+1 : j*M1_cols,:); end
% figure 2
for j = 1 : M2_cols   U3(:,:,j) = OUTPUT_FFTW_REAL((j-1)*M2_cols+1 : j*M2_cols,:); end
for j = 1 : M2_cols   U4(:,:,j) = OUTPUT_FFTW_IMAG((j-1)*M2_cols+1 : j*M2_cols,:); end
for j = 1 : M3_cols   U5(:,:,j) = OUTPUT_INV_FFTW_REAL((j-1)*M3_cols+1 : j*M3_cols,:); end
for j = 1 : M3_cols   U6(:,:,j) = OUTPUT_INV_FFTW_IMAG((j-1)*M3_cols+1 : j*M3_cols,:); end
for j = 1 : M4_cols   U7(:,:,j) = OUTPUT_SWAP_REAL((j-1)*M4_cols+1 : j*M4_cols,:); end
for j = 1 : M4_cols   U8(:,:,j) = OUTPUT_SWAP_IMAG((j-1)*M4_cols+1 : j*M4_cols,:); end
% figure 3
for j = 1 : M5_cols   U9(:,:,j)  = OUTPUT_TRUNC_FFTW_REAL((j-1)*M5_cols+1 : j*M5_cols,:); end
for j = 1 : M5_cols   U10(:,:,j) = OUTPUT_TRUNC_FFTW_IMAG((j-1)*M5_cols+1 : j*M5_cols,:); end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subplot Array Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0;
b=1;

% OPEN FIGURE WINDOW
figure(1);
% % OPEN FIGURE WINDOW
%figure('Color','w');

subplot(1,3,1)
Nx = M3_cols;
Ny = Nx;
Nz = Nx;
s = 0.8; % control little cubes size
ff = .4; % control big cube size
SVL_boxes(ff,a,b,s,Sx,Sy,Sz,Nx,Ny,Nz,U3)
%saveas(gcf,'3D_FFTW.eps')

subplot(1,3,2)
Nx = M4_cols;
Ny = Nx;
Nz = Nx;
s = 0.8;
ff = .4;
SVL_boxes(ff,a,b,s,Sx,Sy,Sz,Nx,Ny,Nz,U7)
% %saveas(gcf,'3D_SWAP.eps')
% 
subplot(1,3,3)
Nx = M5_cols;
Ny = Nx;
Nz = Nx;
s = 0.8;
ff = .43;
SVL_boxes(ff,a,b,s,Sx,Sy,Sz,Nx,Ny,Nz,U9)
% %saveas(gcf,'3D_TRUNCATE_FFTW.eps')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     This is on SVL_boxes.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step size
% dx = Sx/(M3_cols-1); % x-axis step size
% dy = Sy/(M3_cols-1); % y-axis step size
% dz = Sz/(M3_cols-1); % z-axis step size
% 
% xa = 0 : dx : 1;
% ya = 0 : dy : 1;
% za = 0 : dz : 1;
% 
% [Y,X,Z] = meshgrid(ya,xa,za);
% 
% % CREATE RSQ
% %A = X.^2 + Y.^2 + Z.^2;
% A = U3;
% 
% % MASK OF CUBES TO DRAW (THE CUTAWAY)
% MAP = (X - Sx/2).^2 + (Y - Sy/2).^2 + (Z - Sz/2).^2;
% MAP = MAP > (Sx/2)^2;
% 
% % DRAW THE VOXEL VIEW
% %subplot(131);
% hold on;
% 
% % GET METRICS
% Amin = min(U3(:));
% Amax = max(U3(:));
% 
% % SET COLORMAP
% NC   = 256;
% CMAP = parula(NC);
% 
% % DRAW
% s = 0.8;
% 
% for nz = 1 : Nz
%     z1 = za(nz) - s*dz/2;
%     z2 = za(nz) + s*dz/2;
%     
%     for ny = 1 : Ny
%         y1 = ya(ny) - s*dy/2;
%         y2 = ya(ny) + s*dy/2;
%         
%         for nx = 1 : Nx
%             x1 = xa(nx) - s*dx/2;
%             x2 = xa(nx) + s*dx/2;
%             
%             if MAP(nx,ny,nz)
%                 % Color
%                 n = (A(nx,ny,nz) - Amin)/(Amax - Amin);
%                 n = 1 + floor(0.9999*n*NC);
%                 c = CMAP(n,:);
%                 
%                 %xlo
%                 x = x1 * [1 1 1 1 1];
%                 y = [ y1 y2 y2 y1 y1 ];
%                 z = [ z1 z1 z2 z2 z1 ];
%                 fill3(x,y,z,c);
% 
%                 %xhi
%                 x = x2 * [1 1 1 1 1];
%                 y = [ y1 y2 y2 y1 y1 ];
%                 z = [ z1 z1 z2 z2 z1 ];
%                 fill3(x,y,z,c);
% 
%                 %ylo
%                 x = [ x1 x2 x2 x1 x1 ];
%                 y = y1 * [1 1 1 1 1];
%                 z = [ z1 z1 z2 z2 z1 ];
%                 fill3(x,y,z,c);
% 
%                 %yhi
%                 x = [ x1 x2 x2 x1 x1 ];
%                 y = y2 * [1 1 1 1 1];
%                 z = [ z1 z1 z2 z2 z1 ];
%                 fill3(x,y,z,c);
% 
%                 %zlo
%                 x = [ x1 x2 x2 x1 x1 ];
%                 y = [ y1 y1 y2 y2 y1 ];
%                 z = z1 * [1 1 1 1 1];
%                 fill3(x,y,z,c);
% 
%                 %zhi
%                 x = [ x1 x2 x2 x1 x1 ];
%                 y = [ y1 y1 y2 y2 y1 ];
%                 z = z2 * [1 1 1 1 1];
%                 fill3(x,y,z,c);
%             end
%         end
%     end
% end
% 
% % SET GRAPHICS VIEW
% hold off;
% axis equal tight;
% view(-30,30);
% xlabel('$x$','Interpreter','LaTex');
% ylabel('$y$','Interpreter','LaTex');
% zlabel('$z$','Interpreter','LaTex','Rotation',0);
% title('FFTW Voxel Visualization');
% set(gca,'FontSize',18);
