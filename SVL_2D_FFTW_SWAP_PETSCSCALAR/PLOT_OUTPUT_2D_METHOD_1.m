%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables

% Grid size
Lx = 1; % x-axis unit cell grid size
Ly = 1; % y-axis unit cell grid size

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Loading Real Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('OUTPUT_UNIT_CELL_REAL')
%load ('OUTPUT_UNIT_CELL_IMAG')
load ('OUTPUT_FFTW_REAL_X')
load ('OUTPUT_FFTW_IMAG_X')
load ('OUTPUT_FFTW_REAL_Y')
load ('OUTPUT_FFTW_IMAG_Y')
load ('OUTPUT_INV_FFTW_REAL')
load ('OUTPUT_INV_FFTW_IMAG')
load ('OUTPUT_SWAP_REAL_X')
load ('OUTPUT_SWAP_IMAG_X')
load ('OUTPUT_SWAP_REAL_Y')
load ('OUTPUT_SWAP_IMAG_Y')
load ('OUTPUT_TRUNC_FFTW_REAL')
load ('OUTPUT_TRUNC_FFTW_IMAG')
load ('OUTPUT_UNSWAP_FFTW_REAL')
load ('OUTPUT_UNSWAP_FFTW_IMAG')
load ('OUTPUT_INV_FFTW_TRUNC_REAL')
load ('OUTPUT_INV_FFTW_TRUNC_IMAG')

M1 = length(OUTPUT_UNIT_CELL_REAL(:,1));   %  Unit Cell matrix    	(Nx x Ny)
%M2 = length(OUTPUT_UNIT_CELL_IMAG(:,1));   %  Unit Cell matrix    	(Nx x Ny)
M3 = length(OUTPUT_FFTW_REAL_X(:,1)); 	   %  FFTW REAL matrix    	(Nx x Ny)
M4 = length(OUTPUT_FFTW_IMAG_X(:,1)); 	   %  FFTW IMAG matrix    	(Nx x Ny)
M5 = length(OUTPUT_FFTW_REAL_Y(:,1)); 	   %  FFTW REAL matrix    	(Nx x Ny)
M6 = length(OUTPUT_FFTW_IMAG_Y(:,1)); 	   %  FFTW IMAG matrix    	(Nx x Ny)
M7 = length(OUTPUT_INV_FFTW_REAL(:,1));    %  FFTW REAL matrix    	(Nx x Ny)
M8 = length(OUTPUT_INV_FFTW_IMAG(:,1));    %  FFTW IMAG matrix    	(Nx x Ny)
M9 = length(OUTPUT_SWAP_REAL_X(:,1)); 	   %  FFTW REAL matrix    	(Nx x Ny)
M10 = length(OUTPUT_SWAP_IMAG_X(:,1)); 	   %  FFTW IMAG matrix    	(Nx x Ny)
M11 = length(OUTPUT_SWAP_REAL_Y(:,1)); 	   %  FFTW REAL matrix    	(Nx x Ny)
M12 = length(OUTPUT_SWAP_IMAG_Y(:,1)); 	   %  FFTW IMAG matrix    	(Nx x Ny)
M13 = length(OUTPUT_TRUNC_FFTW_REAL(:,1)); %  TRUC FFTW REAL matrix   (NM x NN)
M14 = length(OUTPUT_TRUNC_FFTW_IMAG(:,1)); %  TRUC FFTW IMAG matrix   (NM x NN)
M15 = length(OUTPUT_UNSWAP_FFTW_REAL(:,1));  %  SWAP FFTW REAL matrix   (Nx x Ny)
M16 = length(OUTPUT_UNSWAP_FFTW_IMAG(:,1));  %  SWAP FFTW IMAG matrix   (Nx x Ny)
M17 = length(OUTPUT_INV_FFTW_TRUNC_REAL(:,1));   %  FFTW REAL matrix    	(Nx x Ny)
M18 = length(OUTPUT_INV_FFTW_TRUNC_IMAG(:,1));   %  FFTW IMAG matrix    	(Nx x Ny)

for j = 1 : M1     U1(:,j) =  OUTPUT_UNIT_CELL_REAL(:,j);  end
for j = 1 : M2     U2(:,j) =  OUTPUT_UNIT_CELL_IMAG(:,j);  end
for j = 1 : M3     U3(:,j) =  OUTPUT_FFTW_REAL_X(:,j);     end
for j = 1 : M4     U4(:,j) =  OUTPUT_FFTW_IMAG_X(:,j);     end
for j = 1 : M5     U5(:,j) =  OUTPUT_FFTW_REAL_Y(:,j);     end
for j = 1 : M6     U6(:,j) =  OUTPUT_FFTW_IMAG_Y(:,j);     end
for j = 1 : M7     U7(:,j) =  OUTPUT_INV_FFTW_REAL(:,j);   end
for j = 1 : M8     U8(:,j) =  OUTPUT_INV_FFTW_IMAG(:,j);   end
for j = 1 : M9     U9(:,j) =  OUTPUT_SWAP_REAL_X(:,j);     end
for j = 1 : M10   U10(:,j) =  OUTPUT_SWAP_IMAG_X(:,j);     end
for j = 1 : M11   U11(:,j) =  OUTPUT_SWAP_REAL_Y(:,j);     end
for j = 1 : M12   U12(:,j) =  OUTPUT_SWAP_IMAG_Y(:,j);     end
for j = 1 : M13   U13(:,j) =  OUTPUT_TRUNC_FFTW_REAL(:,j); end
for j = 1 : M14   U14(:,j) =  OUTPUT_TRUNC_FFTW_IMAG(:,j); end
for j = 1 : M15  U15(:,j) =  OUTPUT_UNSWAP_FFTW_REAL(:,j);  end
for j = 1 : M16  U16(:,j) =  OUTPUT_UNSWAP_FFTW_IMAG(:,j);  end
for j = 1 : M17  U17(:,j) =  OUTPUT_INV_FFTW_TRUNC_REAL(:,j);   end
for j = 1 : M18  U18(:,j) =  OUTPUT_INV_FFTW_TRUNC_IMAG(:,j);   end
figure(1)
% Step size
dx = Lx/(M1-1); % x-axis step size
dy = Ly/(M1-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

U1 = rot90(U1,-1); % rotation 90 counterclockwise (-90)
subplot(2,1,1);
imagesc(x, y, U1)      % Unit Cell matrix (Nx x Ny)
axis equal tight
grid on

xlabel('x')
ylabel('y')
title('UNIT CELL REAL')

U2 = rot90(U2,-1); % rotation 90 counterclockwise (-90)
subplot(2,1,2);
imagesc(x, y, U2)      % Unit Cell matrix (Nx x Ny)
axis equal tight
grid on

xlabel('x')
ylabel('y')
title('UNIT CELL IMAG')

% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
print -deps -color OUTPUT_1.eps

figure(2)
% Step size
dx = Lx/(M3-1); % x-axis step size
dy = Ly/(M3-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;


subplot(4,2,1);
imagesc(x, y, U3)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW REAL X')

subplot(4,2,2);
imagesc(x, y, U4)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW IMAG X')

subplot(4,2,3);
imagesc(x, y, U5)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW REAL Y')

subplot(4,2,4);
imagesc(x, y, U6)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW IMAG Y')

subplot(4,2,3);
%U7 = rot90(U7,-1); % rotation 90 counterclockwise (-90)
imagesc(x, y, U7)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('INV FFTW REAL')

subplot(4,2,4);
%U8 = rot90(U8,-1); % rotation 90 counterclockwise (-90)
imagesc(x, y, U8)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('INV FFTW IMAG')

subplot(4,2,5);
imagesc(x, y, U9)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('SWAP REAL X')

subplot(4,2,6);
imagesc(x, y, U10)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('SWAP IMAG X')

subplot(4,2,7);
imagesc(x, y, U11)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('SWAP REAL Y')

subplot(4,2,8);
imagesc(x, y, U12)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('SWAP IMAG Y')

% SAVE PLOTS: saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
print -deps -color OUTPUT_2.eps
%pause

figure(3)
% Step size
dx = Lx/(M13-1); % x-axis step size
dy = Ly/(M13-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

subplot(3,2,1);
imagesc(x, y, U13)      % TRUNC REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW REAL TRUNC')

subplot(3,2,2);
imagesc(x, y, U14)      % TRUNC IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW IMAG TRUNC')
subplot(3,2,3);
imagesc(x, y, U15)      % SWAP REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('UNSWAP TRUNC FFTW REAL')

subplot(3,2,4);
imagesc(x, y, U16)      % SWAP IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('UNSWAP TRUNC FFTW IMAG')

subplot(3,2,5);
%U13 = rot90(U13,-1); % rotation 90 counterclockwise (-90)
imagesc(x, y, U17)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('INV FFTW TRUNC REAL')

subplot(3,2,6);
%U14 = rot90(U14,-1); % rotation 90 counterclockwise (-90)
imagesc(x, y, U18)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('INV FFTW TRUNC IMAG')

% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
print -deps -color OUTPUT_3.eps
disp ("wait please..., hit enter to contnue");
pause (5);

