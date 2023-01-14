function [] = SVL_boxes(ff,a,b,s,Sx,Sy,Sz,Nx,Ny,Nz,A)

% % % OPEN FIGURE WINDOW
% figure('Color','w');

% DRAW THE VOXEL VIEW
hold on;

% ff=.5;
% a=0;
% b=1;

% % Step size
% dx = Sx/(Nx-1); % x-axis step size
% dy = Sy/(Ny-1); % y-axis step size
% dz = Sz/(Nz-1); % z-axis step size
% 
% xa = ff * ( a : dx : b);
% ya = ff * ( a : dx : b);
% za = ff * ( a : dx : b);

xa = ff*Sx*linspace(a,b,Nx);
ya = ff*Sy*linspace(a,b,Ny);
za = ff*Sz*linspace(a,b,Nz);

% increment
dx = xa(2) - xa(1);
dy = ya(2) - ya(1);
dz = za(2) - za(1);

[Y,X,Z] = meshgrid(ya,xa,za);

% MASK OF CUBES TO DRAW (THE CUTAWAY)
MAP = (X - Sx/2).^2 + (Y - Sy/2).^2 + (Z - Sz/2).^2;
MAP = MAP > (Sx/2)^2;

% GET METRICS
Amin = min(A(:));
Amax = max(A(:));

% SET COLORMAP
NC   = 256;
CMAP = parula(NC);

for nz = 1 : Nz
    z1 = za(nz) - s*dz/2;
    z2 = za(nz) + s*dz/2;
    
    for ny = 1 : Ny
        y1 = ya(ny) - s*dy/2;
        y2 = ya(ny) + s*dy/2;
        
        for nx = 1 : Nx
            x1 = xa(nx) - s*dx/2;
            x2 = xa(nx) + s*dx/2;
            
            if MAP(nx,ny,nz)
                % Color
                n = (A(nx,ny,nz) - Amin)/(Amax - Amin);
                n = 1 + floor(0.9999*n*NC);
                c = CMAP(n,:);
                
                %xlo
                x = x1 * [1 1 1 1 1];
                y = [ y1 y2 y2 y1 y1 ];
                z = [ z1 z1 z2 z2 z1 ];
                fill3(x,y,z,c);

                %xhi
                x = x2 * [1 1 1 1 1];
                y = [ y1 y2 y2 y1 y1 ];
                z = [ z1 z1 z2 z2 z1 ];
                fill3(x,y,z,c);

                %ylo
                x = [ x1 x2 x2 x1 x1 ];
                y = y1 * [1 1 1 1 1];
                z = [ z1 z1 z2 z2 z1 ];
                fill3(x,y,z,c);

                %yhi
                x = [ x1 x2 x2 x1 x1 ];
                y = y2 * [1 1 1 1 1];
                z = [ z1 z1 z2 z2 z1 ];
                fill3(x,y,z,c);

                %zlo
                x = [ x1 x2 x2 x1 x1 ];
                y = [ y1 y1 y2 y2 y1 ];
                z = z1 * [1 1 1 1 1];
                fill3(x,y,z,c);

                %zhi
                x = [ x1 x2 x2 x1 x1 ];
                y = [ y1 y1 y2 y2 y1 ];
                z = z2 * [1 1 1 1 1];
                fill3(x,y,z,c);
            end
        end
    end
end

% SET GRAPHICS VIEW
%hold off;
axis equal tight;
view(120,30);
xlabel('$p$','Interpreter','LaTex');
ylabel('$q$','Interpreter','LaTex');
zlabel('$r$','Interpreter','LaTex','Rotation',0);
%title('Voxel Visualization');
set(gca,'FontSize',25);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
%colormap hot %(parula(5))%winter
colorbar
%set(gca,'xticklabel',0:0.2:1)
%axis([a b a b])


