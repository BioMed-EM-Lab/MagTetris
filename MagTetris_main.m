% Examples of MagTetris Calculation
% Each section is an example of MagTetris function. They can be executed individually.
% @author  Ting-Ou Liang
% @version 2024/03/11

%% Section 1: Field calculation and visualization of cuboid magnets
clear;
close all;

% Field-of-view in [mm]
FOV_half_1 = 50;                   % Half of the length for the first dimension
FOV_1 = -FOV_half_1:1:FOV_half_1;

FOV_half_2 = 50;                   % Half of the length for the second dimension
FOV_2 = -FOV_half_2:1:FOV_half_2;

FOV_3 = 0;                          % The location of the plane in the third dimension

surface = 'z';      % 'x' - yz plane, 'y' - xz plane, 'z' - xy plane


% The definition of the three cuboid magnets
% Center location (x,y,z) [mm]
loc_all_list = [80,0,0;0,50,-50;-45,-45,0];
% Orientation angle (yaw,pitch,roll) [degree]
angle_all = [0,0,0;45,0,0;0,0,90];
% Magnet dimensions (x,y,z) [mm]
magnet_dim_all = [50,20,30;50,50,50;20,20,20];
% Remanence in [T]
Br_all = [1.40,1.43,1.38];                   

Three_MT = MagTetris();
Three_MT = Three_MT.AssignCuboid(loc_all_list,angle_all,Br_all,magnet_dim_all);
n_per_group = 1;    % The number of magnets for grouping during calculation, by default it should be 1
[Bx,By,Bz] = Three_MT.Field2D(FOV_1,FOV_2,FOV_3,surface,n_per_group);
% Only keep the circular region
for idx_1=1:length(FOV_1)
    for idx_2=1:length(FOV_2)
        r = sqrt(FOV_1(idx_1)^2 + FOV_2(idx_2)^2);
        if r > FOV_half_1
            Bx(idx_2,idx_1) = NaN;
            By(idx_2,idx_1) = NaN;
            Bz(idx_2,idx_1) = NaN;
        end
    end
end

% Change the unit from [T] to [mT]
Bx = Bx*1e3;
By = By*1e3;
Bz = Bz*1e3;

% Plot the Bx/By/Bz field component
font_size = 18;
figure(1);
pcolor(FOV_1,FOV_2,Bx);
axis square;
shading flat;
cb = colorbar;
set(get(cb,'Title'),'string','mT');
colormap jet;
set(gcf,'color','w');
ax = gca;
ax.FontSize = font_size;
if surface == 'x'
    xlabel('y/mm','FontSize',font_size)
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('Bx at x = %.1f mm (mT)',FOV_3);
elseif surface == 'y'
    xlabel('x/mm','FontSize',font_size)
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('Bx at y = %.1f mm (mT)',FOV_3);
elseif surface == 'z'
    xlabel('x/mm','FontSize',font_size);
    ylabel('y/mm','FontSize',font_size);
    ttl = sprintf('Bx at z = %.1f mm (mT)',FOV_3);
end
title(ttl);

figure(2);
pcolor(FOV_1,FOV_2,By);
axis square;
shading flat;
cb = colorbar;
set(get(cb,'Title'),'string','mT');
colormap jet;
set(gcf,'color','w');
ax = gca;
ax.FontSize = font_size;
if surface == 'x'
    xlabel('y/mm','FontSize',font_size)
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('By at x = %.1f mm (mT)',FOV_3);
elseif surface == 'y'
    xlabel('x/mm','FontSize',font_size)
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('By at y = %.1f mm (mT)',FOV_3);
elseif surface == 'z'
    xlabel('x/mm','FontSize',font_size);
    ylabel('y/mm','FontSize',font_size);
    ttl = sprintf('By at z = %.1f mm (mT)',FOV_3);
end
title(ttl);

figure(3);
pcolor(FOV_1,FOV_2,Bz);
axis square;
shading flat;
cb = colorbar;
set(get(cb,'Title'),'string','mT');
colormap jet;
set(gcf,'color','w');
ax = gca;
ax.FontSize = font_size;
if surface == 'x'
    xlabel('y/mm','FontSize',font_size)
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('Bz at x = %.1f mm (mT)',FOV_3);
elseif surface == 'y'
    xlabel('x/mm','FontSize',font_size)
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('Bz at y = %.1f mm (mT)',FOV_3);
elseif surface == 'z'
    xlabel('x/mm','FontSize',font_size);
    ylabel('y/mm','FontSize',font_size);
    ttl = sprintf('Bz at z = %.1f mm (mT)',FOV_3);
end
title(ttl);

% Visualization of the magnetization
DrawMagnetization(loc_all_list,angle_all,magnet_dim_all,FOV_1,FOV_2,FOV_3,surface,Bz);

% Weight of the PMA
weight = Three_MT.Weight();
fprintf('Weight: %.1f kg\n', weight);

%% Section 2: Force calculation of 2 cuboid magnets
% Please see the ForceSingle function in MagTetris class for details
% The computed forces have the unit [N]

% The definition of the two cuboid magnets
loc_all = [0,0,0;0,25,-15];
angle_all = [0,0,0;180,0,0];
dim_all = [10,10,10;10,10,10];
Br_all = [1.40,1.43];

two_cub_MT = MagTetris();
% Cuboid
two_cub_MT = two_cub_MT.AssignCuboid(loc_all,angle_all,Br_all,dim_all);

% The number of discretization segments for each dimension
% The surface of magnet will be divided into small rectangles for force
% calculation. Larger FOV_del_ratio will result in more accurate estimation.
FOV_del_ratio = 20;

% Magnetic force acting on magnet #1 from magnet #2
% The [1,1] is the target, and the format is [(type of magnet),(index of magnet)].
% Type = 1 refers to cuboid magnet. Currently the code only support cuboid magnets.
F_single_1 = two_cub_MT.ForceSingle([1,1],[1,2],FOV_del_ratio);
fprintf('Force from #2 to #1: (Fx,Fy,Fz) = (%.1f,%.1f,%.1f) N\n', F_single_1);
% Magnetic force acting on magnet #2 from magnet #1
F_single_2 = two_cub_MT.ForceSingle([1,2],[1,1],FOV_del_ratio);
fprintf('Force from #1 to #2: (Fx,Fy,Fz) = (%.1f,%.1f,%.1f) N\n', F_single_1);
% The forces experienced by all magnets within the PMA
F_all_static = two_cub_MT.ForceAllStatic(FOV_del_ratio);
for magnet_idx=1:size(F_all_static,1)
    fprintf('Force acting on magnet #%d: (Fx,Fy,Fz) = (%.1f,%.1f,%.1f) N\n', magnet_idx,F_all_static(magnet_idx,:));
end
