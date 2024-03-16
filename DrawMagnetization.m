% This function helps to visualize cuboid magnets and their magnetization.
% @author  Junqi Yang, Ting-Ou Liang
% @version 2024/02/28


function DrawMagnetization(loc_all,angle_all,magnet_dim_all,FOV_1,FOV_2,FOV_3,surface,B_field)
% INPUT:
%       loc_all_list - the center locations of all cuboid
%       angle_all - the orientation angle (yaw,pitch,roll) with respect to
%                   the center of each magnet
%       magnet_dim_all - the (len_x,len_y,len_z) dimension of the magnets
%                        without rotation
%       FOV_1/2[1,m1]/[1,m2] - FOV 1D arrays of two directions
%       FOV_3 - a scalar, the third dimension of field calculation
%       surface - z: XY plane, y: XZ plane, x: YZ plane, where ab-plane means FOV_1 = FOV_a, FOV_2 = FOV_b
%       B_field - the 2D magnetic field produced by the PMA

figure;
view(3); 
hold on;
axis equal;
grid on;
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');
title('3D Visualization of Magnets');
set(gcf,'color','w');

for i = 1:size(loc_all, 1)
    center = loc_all(i, :);
    angles = angle_all(i, :); % Yaw, Pitch, Roll in degrees
    dimensions = magnet_dim_all(i, :); % Dimensions [x, y, z]

    DrawRotatedCuboid(center,angles,dimensions);
end
hold on;

% Vector initialization
dir_vectors = zeros(size(loc_all));
for i = 1:size(loc_all, 1)
    % Create rotation matrix
    R_yaw = [cosd(angle_all(i,1)) -sind(angle_all(i,1)) 0; sind(angle_all(i,1)) cosd(angle_all(i,1)) 0; 0 0 1];
    R_pitch = [cosd(angle_all(i,2)) 0 sind(angle_all(i,2)); 0 1 0; -sind(angle_all(i,2)) 0 cosd(angle_all(i,2))];
    R_roll = [1 0 0; 0 cosd(angle_all(i,3)) -sind(angle_all(i,3)); 0 sind(angle_all(i,3)) cosd(angle_all(i,3))];
    
    % Perform the rotation
    dir_vector = R_yaw * R_pitch * R_roll * [0; 1; 0];
    dir_vectors(i, :) = dir_vector';
end

% Rescale the length of the vector
% Decide the scale factor by the dimension of magnets
scale_factor = 2*magnet_dim_all(:,2);
scaled_dir_vectors = dir_vectors .* scale_factor;

% Compute the new starting point of the vectors
new_start_points = loc_all - 0.5 * scaled_dir_vectors;

% Draw the vectors with adjusted locations and orientations
quiver3(new_start_points(:,1), new_start_points(:,2), new_start_points(:,3), scaled_dir_vectors(:,1), scaled_dir_vectors(:,2), scaled_dir_vectors(:,3), 'AutoScale', 'off');
hold on;

% Visualize the magnetic field on the observation plane
switch surface
    case 'x'
        [Y, Z] = meshgrid(linspace(min(FOV_1(:)),max(FOV_1(:)),size(B_field, 2)),linspace(min(FOV_2(:)),max(FOV_2(:)),size(B_field, 1)));
        X = FOV_3*ones(size(Y)); 
    case 'y'
        [X, Z] = meshgrid(linspace(min(FOV_1(:)),max(FOV_1(:)),size(B_field, 2)),linspace(min(FOV_2(:)),max(FOV_2(:)),size(B_field, 1)));
        Y = FOV_3*ones(size(X)); 
    case 'z'
        [X, Y] = meshgrid(linspace(min(FOV_1(:)),max(FOV_1(:)),size(B_field, 2)),linspace(min(FOV_2(:)),max(FOV_2(:)),size(B_field, 1)));
        Z = FOV_3*ones(size(X)); 
end

surf(X,Y,Z,B_field,'EdgeColor','none');
colormap jet; 
cb = colorbar;
set(get(cb,'Title'),'string','mT');
end

%% A helper function to draw cuboids
function DrawRotatedCuboid(center,angles,dimensions)
% Generate the vertices
[X, Y, Z] = meshgrid([-1, 1], [-1, 1], [-1, 1]);
vertices = [X(:), Y(:), Z(:)] .* 0.5;
vertices = vertices .* dimensions;
vertices = vertices + center;

% Apply the rotation
R_yaw = [cosd(angles(1)),-sind(angles(1)),0;...
         sind(angles(1)),cosd(angles(1)),0;...
         0,0,1];
R_pitch = [cosd(angles(2)),0,sind(angles(2));...
           0,1,0;...
           -sind(angles(2)),0,cosd(angles(2))];
R_roll = [1,0,0;...
          0,cosd(angles(3)),-sind(angles(3));...
          0,sind(angles(3)),cosd(angles(3))];
R = R_yaw*R_pitch*R_roll;

% Shift and rotate the vertices
vertices = (R * (vertices - center)')' + center;
% Define the faces of the cuboid
faces = [1,2,4,3; 5,6,8,7; 1,5,7,3; 2,6,8,4; 1,5,6,2; 3,7,8,4];
% Draw the cuboids
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'red', 'EdgeColor', 'k', 'LineWidth', 1.5);

end
