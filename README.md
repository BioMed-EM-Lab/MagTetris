
# MagTetris: A simulator for fast magnetic field and force calculation for permanent magnet array designs

`MagTetris` is a simulator that
can perform fast calculations for both the B-field and the magnetic
force simultaneously for a PMA consisting of cuboid magnets.

![Logo](https://github.com/TingouLiang/MagTetris/blob/main/MagTetris%20Logo.png?raw=true)


## Related Publication
You can check the following publication to learn the framework and performance of MagTetris:

Ting-Ou Liang, Yan Hao Koh, Tie Qiu, Erping Li, Wenwei Yu, Shao Ying Huang,
*MagTetris: A simulator for fast magnetic field and force calculation for permanent magnet array designs,
Journal of Magnetic Resonance*,
Volume 352,
2023,
107463,
ISSN 1090-7807,
https://doi.org/10.1016/j.jmr.2023.107463.
## File List
- `MagTetris.m` - the complete class definition of `MagTetris`
- `MagTetris_main.m` - the scripts that contains examples for each `MagTetris` feature
- `DrawMagnetization.m` - a function to visualize the magnets with their magnetization and the magnetic field on the observation plane

For the usage of magnetic field/force calculation and other functions, please check `MagTetris_main.m` for the details. Each section is an individual demo for different features.
## Usage/Examples
Here is a section taken out from `MagTetris_main.m`, which shows the basic functions of magnetic field calculation for cuboid magnets, etc.

An example of force calculation is in Section 2 of the script `MagTetris_main.m`. Please find more detail in the file. The syntax and usage are explained in the comments of the code.

A step-to-step explanation of the first section is shown below:

### Step 1: Field-of-view definition
The observation plane is described by the following quantities:
- `FOV_1` - the horizontal dimension of the observation plane
- `FOV_2` - the vertical dimension of the observation plane
- `FOV_3` - a single value to indicate the location of the plane along the third dimension
- `surface` - the orientation of the observation plane: `'x'` - yz plane, `'y'` - xz plane, `'z'` - xy plane

```matlab
%% Field calculation: cuboid magnets
clear;
close all;

% Field-of-view in [mm]
FOV_half_1 = 50;                   % Half of the length for the first dimension
FOV_1 = -FOV_half_1:1:FOV_half_1;

FOV_half_2 = 50;                   % Half of the length for the second dimension
FOV_2 = -FOV_half_2:1:FOV_half_2;

FOV_3 = 0;                          % The location of the plane in the third dimension

surface = 'z';      % 'x' - yz plane, 'y' - xz plane, 'z' - xy plane
```

### Step 2: Permanet magnet array definition
A permanent magnet array consisting of cuboid magnets is defined by four parameters:
- `loc_all_list` - the center location of each magnet in the form of `(x, y, z)` in row-order
- `angle_all` - the `(yaw,pitch,roll)` orientation angle with respect to the center of each magnet. An angle of `(0,0,0)` indicates +y-magnetization, and cuboid sides should be parallel with the Cartesian coordinates.
- `magnet_dim_all` - the side length of the cuboid `(len_x,len_y,len_z)` for all the magnets. The side length is measured when the magnet is placed with an orientation angle of `(0,0,0)`
- `Br_all` - the magnetic remanence

In this example, three cuboid magnets are created.

```matlab
% The definition of the three cuboid magnets
% Center location (x,y,z) [mm]
loc_all_list = [80,0,0;0,50,-50;-45,-45,0];
% Orientation angle (yaw,pitch,roll) [degree]
angle_all = [0,0,0;45,0,0;0,0,90];
% Magnet dimensions (x,y,z) [mm]
magnet_dim_all = [50,50,20;50,50,50;20,20,20];
% Remanence in [T]
Br_all = [1.40,1.43,1.38];                   
```

### Step 3: Magnetic field calculation
Create a MagTetris object using the function `MagTetris()`, and use the method `AssignCuboid` to input the three-magnet array defined above.
The method `Field2D` is used to compute the `Bx`, `By`, and `Bz` component of the magnetic field within the observation plane.

```matlab
Three_MT = MagTetris();
Three_MT = Three_MT.AssignCuboid(loc_all_list,angle_all,Br_all,magnet_dim_all);
n_per_group = 1;    % The number of magnets for grouping during calculation, by default it should be 1
[Bx,By,Bz] = Three_MT.Field2D(FOV_1,FOV_2,FOV_3,surface,n_per_group);
```

### Step 4: Post-processing of the result
In this example, the observation region is circular, and the displayed unit is changed from [T] to [mT].

```matlab
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
```

### Step 5: Plot the 2D magnetic field maps
The field components in x-, y-, and z-direction are plotted in separate figures. The output will be three figures as shown below:
![Field](https://github.com/TingouLiang/MagTetris/blob/main/MagTetris%20Field.png?raw=true)

```matlab
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
xlabel('x/mm','FontSize',font_size);
if surface == 'y'
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('Bx at y = %.1f mm (mT)',FOV_3);
else
    ylabel('y/mm','FontSize',font_size);
    ttl = sprintf('Bx at z = %.1f mm (mT)',FOV_3);
end

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
xlabel('x/mm','FontSize',font_size);
if surface == 'y'
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('By at y = %.1f mm (mT)',FOV_3);
else
    ylabel('y/mm','FontSize',font_size);
    ttl = sprintf('By at z = %.1f mm (mT)',FOV_3);
end

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
xlabel('x/mm','FontSize',font_size);
if surface == 'y'
    ylabel('z/mm','FontSize',font_size);
    ttl = sprintf('Bz at y = %.1f mm (mT)',FOV_3);
else
    ylabel('y/mm','FontSize',font_size);
    ttl = sprintf('Bz at z = %.1f mm (mT)',FOV_3);
end
```

### Step 6: Visualization of the magnetization
A visualization of the magnets can be performed using the function `DrawMagnetization`. The corresponding magnetizations will be drawn accordingly as well. The calculated magnetic field can be shown on the observation plane. The output will be a figure as shown below:
![Visualization](https://github.com/TingouLiang/MagTetris/blob/main/MagTetris%20Visualization.png?raw=true)

```matlab
% Visualization of the magnetization
DrawMagnetization(loc_all_list,angle_all,magnet_dim_all,FOV_1,FOV_2,FOV_3,surface,Bz);
```

### Step 7: A quick estimation of the weight of the PMA
Finally, `MagTetris` can provide a quick estimation of the weight for the magnet array based on the total volume of the magnets. The reference is the NdFeB magnet, which has a weight of 8g with a dimension of 10\*10\*10 mm^3.

```matlab
% Weight of the PMA
weight = Three_MT.Weight();
fprintf('Weight: %.1f kg\n', weight);

```

