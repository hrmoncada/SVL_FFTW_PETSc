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
load ('OUTPUT_FFTW_REAL')
load ('OUTPUT_FFTW_IMAG')
load ('OUTPUT_INV_FFTW_REAL')
load ('OUTPUT_INV_FFTW_IMAG')
load ('OUTPUT_SWAP_REAL')
load ('OUTPUT_SWAP_IMAG')
load ('OUTPUT_TRUNC_FFTW_REAL')
load ('OUTPUT_TRUNC_FFTW_IMAG')

M1  = length(OUTPUT_UNIT_CELL_REAL(:,1));  %  Unit Cell matrix    	(Nx x Ny)
M3  = length(OUTPUT_FFTW_REAL(:,1)); 	   %  FFTW REAL matrix    	(Nx x Ny)
M4  = length(OUTPUT_FFTW_IMAG(:,1)); 	   %  FFTW IMAG matrix    	(Nx x Ny)
M5  = length(OUTPUT_INV_FFTW_REAL(:,1));   %  FFTW REAL matrix    	(Nx x Ny)
M6  = length(OUTPUT_INV_FFTW_IMAG(:,1));   %  FFTW IMAG matrix    	(Nx x Ny)
M7  = length(OUTPUT_SWAP_REAL(:,1));  %  SWAP FFTW REAL matrix   (Nx x Ny)
M8  = length(OUTPUT_SWAP_IMAG(:,1));  %  SWAP FFTW IMAG matrix   (Nx x Ny)
M9  = length(OUTPUT_TRUNC_FFTW_REAL(:,1)); %  TRUC FFTW REAL matrix   (NM x NN)
M10 = length(OUTPUT_TRUNC_FFTW_IMAG(:,1)); %  TRUC FFTW IMAG matrix   (NM x NN)

for j = 1 : M1    U1(:,j) =  OUTPUT_UNIT_CELL_REAL(:,j);  end
for j = 1 : M3    U3(:,j) =  OUTPUT_FFTW_REAL(:,j);       end
for j = 1 : M4    U4(:,j) =  OUTPUT_FFTW_IMAG(:,j);       end
for j = 1 : M5    U5(:,j) =  OUTPUT_INV_FFTW_REAL(:,j);   end
for j = 1 : M6    U6(:,j) =  OUTPUT_INV_FFTW_IMAG(:,j);   end
for j = 1 : M7    U7(:,j) =  OUTPUT_SWAP_REAL(:,j);  end
for j = 1 : M8    U8(:,j) =  OUTPUT_SWAP_IMAG(:,j);  end
for j = 1 : M9    U9(:,j) =  OUTPUT_TRUNC_FFTW_REAL(:,j); end
for j = 1 : M10  U10(:,j) =  OUTPUT_TRUNC_FFTW_IMAG(:,j); end

figure(1)
% Step size
dx = Lx/(M1-1); % x-axis step size
dy = Ly/(M1-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

U1 = rot90(U1,-1); % rotation 90 counterclockwise (-90)
%subplot(1,2,1);
imagesc(x, y, U1)      % Unit Cell matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('UNIT CELL REAL')

% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_1.eps

figure(2)
% Step size
dx = Lx/(M3-1); % x-axis step size
dy = Ly/(M3-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

subplot(1,2,1);
imagesc(x, y, U3)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW REAL')

subplot(1,2,2);
imagesc(x, y, U4)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('FFTW IMAG')

% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_2.eps

figure(3)
% Step size
dx = Lx/(M3-1); % x-axis step size
dy = Ly/(M3-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

subplot(1,2,1);
%U5 = rot90(U5,-1); % rotation 90 counterclockwise (-90)
imagesc(x, y, U5)      % REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('INV FFTW REAL')

subplot(1,2,2);
U6 = rot90(U6,-1); % rotation 90 counterclockwise (-90)
imagesc(x, y, U6)      % IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
xlabel('x')
ylabel('y')
title('INV FFTW IMAG')

% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_3.eps

figure(4)
% Step size
dx = Lx/(M8-1); % x-axis step size
dy = Ly/(M8-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

subplot(2,2,1);
imagesc(x, y, U7)      % SWAP REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
colorbar()
xlabel('x')
ylabel('y')
title('SWAP FFTW REAL')

subplot(2,2,2);
imagesc(x, y, U8)      % SWAP IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
colorbar()
xlabel('x')
ylabel('y')
title('SWAP FFTW IMAG')

dx = Lx/(M9-1); % x-axis step size
dy = Ly/(M9-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;

subplot(2,2,3);
imagesc(x, y, U9)      % TRUNC REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
colorbar()
xlabel('x')
ylabel('y')
title('FFTW REAL TRUNC')

subplot(2,2,4);
imagesc(x, y, U10)      % TRUNC IMAG FFTW matrix (Nx x Ny)
axis equal tight
grid on
colorbar()
xlabel('x')
ylabel('y')
title('FFTW IMAG TRUNC')


% SAVE PLOTS:  saveas (1,"test.eps")  or print (1,"test.eps") or print -deps test.eps
%saveas (2,"test1.eps")
%print (2,"test2.eps")
%print('figure.eps','-deps')
%print -deps -color OUTPUT_4.eps


figure(4)
subplot(1,2,1);
imagesc(x, y, U7)      % SWAP REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
%colorbar()
xlabel('$p$','Interpreter','LaTex');
ylabel('$q$','Interpreter','LaTex','rotation',0);
set(gca,'FontSize',25);
%title('SWAP FFTW REAL')


subplot(1,2,2);
imagesc(x, y, U9)      % TRUNC REAL FFTW matrix (Nx x Ny)
axis equal tight
grid on
colorbar()
xlabel('$p$','Interpreter','LaTex');
ylabel('$q$','Interpreter','LaTex','rotation',0);
set(gca,'FontSize',25);
%title('FFTW REAL TRUNC')

%fprintf(stderr,"\nWait Please...\n");
%disp ("Hit Enter to Continue...\n");
pause;
close all
