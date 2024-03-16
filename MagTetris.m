% The class definition of MagTetris with the support of cuboid permanent magnets.
% @author  Ting-Ou Liang
% @version 2024/03/11


classdef MagTetris

    properties
        % Parameters for cuboid magnets
        cub_loc         % Center location [mm]
        cub_angle       % Orientation angle (yaw,pitch,roll) [degree] with respect to +y polarization
        cub_dim         % Magnet dimension (x,y,z) [mm] when placed to have +y polarization
        cub_Br          % Magnet remanence in [T]

        % Other useful parameters
        % Total number of magnets
        cub_numMag      % The number of cuboid magnets
        
        % Variables for PMA rotation
        cub_locRot      % The center location of all cuboid magnets after a rotation of the whole PMA
        angleRot        % The orientation angle of all magnets after a rotation of the whole PMA
    end

    methods
        % ========= Initialization Functions ========== %
        % Initialization function
        function self = MagTetris()
            self.cub_numMag = 0;
        end

        % Assign design parameters for cuboid magnets
        function self = AssignCuboid(self,loc_list,angle_list,Br_list,dim_list)
        % INPUT:
        %       loc_list - center location of all magnets
        %       angle_list - 3D orientation angles of all magnets
        %       dim_list - dimension of all magnets
        %       Br_list - remanence of all magnets

            % loc: [n,3] matrix with each row as (x,y,z)
		    self.cub_loc = loc_list;
		    % angle: [n,3] matrix with each row as (yaw,pitch,roll)
		    self.cub_angle = angle_list;
		    % dim: [n,3] matrix with each row as (dim_x,dim_y,dim_z)
		    self.cub_dim = dim_list;
		    % Br: length n vector
		    self.cub_Br = Br_list;

            self.cub_numMag = length(self.cub_Br);

            self.cub_locRot = self.cub_loc;
            self.angleRot = self.cub_angle;
        end

        
        % Rotate the whole array w.r.t. a pivot.
        % If not given, the pivot will be the origin (0,0,0)
        % Currently the rotation is only for cuboid magnets
        function self = RotateArray(self,angle_rot,pivot)
        % INPUT:
	    %       angle_rot - (yaw,pitch,roll) for array rotation w.r.t. origin

            % Assign new angle
		    self.angleRot(:,1) = self.cub_angle(:,1) + angle_rot(1);
		    self.angleRot(:,2) = self.cub_angle(:,2) + angle_rot(2);
		    self.angleRot(:,3) = self.cub_angle(:,3) + angle_rot(3);


            % If a pivot is given
            if nargin == 3
                % Shift the location
                loc_shift = self.cub_loc - pivot;
                loc_shift_rot = self.Rotate3D(loc_shift,angle_rot,1);
                self.cub_locRot = loc_shift_rot + pivot;
            else
		        self.cub_locRot = self.Rotate3D(self.cub_loc,angle_rot,1);
            end
        end

        % Concatenate another MagTetris object with the existing one
        function self = ConcatMT(self,other_MT)
        % INPUT:
	    %       other_MT - another MagTetris object to be concatenated
            if other_MT.cub_numMag ~= 0
                self.cub_numMag = self.cub_numMag + other_MT.cub_numMag;
                self.cub_loc = [self.cub_loc;other_MT.cub_loc];
                self.cub_angle = [self.cub_angle;other_MT.cub_angle];
                self.cub_Br = [self.cub_Br,other_MT.cub_Br];
                self.cub_dim = [self.cub_dim;other_MT.cub_dim];
            end

        end
        % ========= Magnetic Field Calculation ========== %

        % Compute magnetic field for a specified 2D FoV
        function [Bx,By,Bz] = Field2D(self,FOV_1,FOV_2,FOV_3,surface,n_per_group)
        % INPUT:
	    %       FOV_1/2[1,m1]/[1,m2] - FOV 1D arrays of two directions
	    %       FOV_3 - a scalar, the third dimension of field calculation
	    %       surface - z: XY plane, y: XZ plane, x: YZ plane, where ab-plane means FOV_1 = FOV_a, FOV_2 = FOV_b
	    %       n_per_group - the number of magnets in a group for multi-magnet method
        %
	    % OUTPUT:
	    %       Bx/By/Bz - the x/y/z field components in (T)

            FOV1_len = length(FOV_1);
            FOV2_len = length(FOV_2);
            Bx = zeros(FOV2_len,FOV1_len);
            By = zeros(FOV2_len,FOV1_len);
            Bz = zeros(FOV2_len,FOV1_len);

            % ===== cuboid magnet =====
            % Group selection
            if n_per_group > 1 && self.cub_numMag ~= 0
                % Multiple magnets in a group
                num_whole = floor(self.cub_numMag/n_per_group);
                n_part = mod(self.cub_numMag,n_per_group);
                for i=0:num_whole
                    idx_start = i*n_per_group + 1;
                    if i == num_whole
                        loc_curr = self.cub_loc(idx_start:end,:);
                        angle_curr = self.cub_angle(idx_start:end,:);
                        dim_curr = self.cub_dim(idx_start:end,:);
                        Br_curr = self.cub_Br(idx_start:end);

                        [Bx_curr, By_curr, Bz_curr] = self.field2DMultiCuboid(FOV1_len,FOV2_len,n_part,FOV_1,FOV_2,FOV_3,dim_curr,loc_curr,angle_curr,Br_curr,surface);

                        Bx = Bx + Bx_curr;
				        By = By + By_curr;
				        Bz = Bz + Bz_curr;

                        %fprintf('Cuboid Magnet #%d\n',self.cub_numMag);
                    else
                        idx_end = (i+1)*n_per_group;
                        loc_curr = self.cub_loc(idx_start:idx_end,:);
                        angle_curr = self.cub_angle(idx_start:idx_end,:);
                        dim_curr = self.cub_dim(idx_start:idx_end,:);
                        Br_curr = self.cub_Br(idx_start:idx_end);

                        [Bx_curr, By_curr, Bz_curr] = self.field2DMultiCuboid(FOV1_len,FOV2_len,n_per_group,FOV_1,FOV_2,FOV_3,dim_curr,loc_curr,angle_curr,Br_curr,surface);

                        Bx = Bx + Bx_curr;
				        By = By + By_curr;
				        Bz = Bz + Bz_curr;

                        %fprintf('Magnet #%d\n',idx_end);
                    end
                end
            else
                % Loop through every magnet
                for mag_idx=1:self.cub_numMag
                    loc_curr = self.cub_loc(mag_idx,:);
                    angle_curr = self.cub_angle(mag_idx,:);
                    dim_curr = self.cub_dim(mag_idx,:);
                    Br_curr = self.cub_Br(mag_idx);

                    [Bx_curr, By_curr, Bz_curr] = self.Field2DSingleCuboid(FOV1_len,FOV2_len,FOV_1,FOV_2,FOV_3,dim_curr,loc_curr,angle_curr,Br_curr,surface);

                    Bx = Bx + Bx_curr;
			        By = By + By_curr;
			        Bz = Bz + Bz_curr;

%                     if mod(mag_idx,50) == 0 || mag_idx == self.cub_numMag
%                         fprintf('Cuboid Magnet #%d\n',mag_idx);
%                     end
                end
            end
        end

        % ========= Helper Functions ========= %

        % Compute the 2D magnetic field produced by a single cuboid magnet
        function [Bx,By,Bz] = Field2DSingleCuboid(self,FOV1_len,FOV2_len,FOV_1,FOV_2,FOV_3,magnet_dim,magnet_loc,angle,Br,surface)
	    % INPUT:
        %       FOV1/2_len - the length of FOV_1/2 vector
	    %       FOV_1/2[1,m1]/[1,m2] - FOV 1D arrays of two directions
        %                        For curvy surface, FOV_1 = phi, FOV_2 = z, FOV_3 = rho
	    %       FOV_3 - a scalar, the third dimension of field calculation
	    %       magnet_dim - the size of the magnet
	    %       magnet_loc - the center location of the magnet
	    %       angle - the orientation angle of magnet in the form of (yaw,pitch,row)
	    %       Br - remanence of magnet
	    %       surface - z: XY plane, y: XZ plane, x: YZ plane, c: curvy,
        %                 where ab-plane means FOV_1 = FOV_a, FOV_2 = FOV_b
        %
	    % OUTPUT:
	    %       Bx/By/Bz - the x/y/z field components in (T)
            
            B_all_local = zeros(FOV1_len*FOV2_len,3);
            if Br == 0
			    Bx = zeros(FOV2_len,FOV1_len);
			    By = zeros(FOV2_len,FOV1_len);
			    Bz = zeros(FOV2_len,FOV1_len);
			    return;
            end

            % Change the coordinates
            % normal vector = y: XZ slice, z: XY slice, x: YZ slice
            if surface == 'x'
                [FOV_y,FOV_z] = meshgrid(FOV_1,FOV_2);
                shifted_x = FOV_3 - magnet_loc(1);
                shifted_y = FOV_y - magnet_loc(2);
                shifted_z = FOV_z - magnet_loc(3);
                FOV_shifted = ones(numel(FOV_y),3);
                FOV_shifted(:,1) = FOV_shifted(:,1)*shifted_x;
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            elseif surface == 'y'
                [FOV_x,FOV_z] = meshgrid(FOV_1,FOV_2);
                shifted_x = FOV_x - magnet_loc(1);
                shifted_y = FOV_3 - magnet_loc(2);
                shifted_z = FOV_z - magnet_loc(3);
                FOV_shifted = ones(numel(FOV_x),3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = FOV_shifted(:,2)*shifted_y;
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            elseif surface == 'z'
                [FOV_x,FOV_y] = meshgrid(FOV_1,FOV_2);
                shifted_x = FOV_x - magnet_loc(1);
                shifted_y = FOV_y - magnet_loc(2);
                shifted_z = FOV_3 - magnet_loc(3);
                FOV_shifted = ones(numel(FOV_x),3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = FOV_shifted(:,3)*shifted_z;
            elseif surface == 'c'
                % Cylindrical surface
                [FOV_phi,FOV_z] = meshgrid(FOV_1,FOV_2);
                FOV_rho = FOV_3*ones(size(FOV_phi));
                [FOV_x,FOV_y] = pol2cart(FOV_phi+pi/2,FOV_rho);
                shifted_x = FOV_x - magnet_loc(1);
                shifted_y = FOV_y - magnet_loc(2);
                shifted_z = FOV_z - magnet_loc(3);
                FOV_shifted = ones(numel(FOV_x),3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            elseif surface == 's'
                % Spherical surface
                [FOV_azi,FOV_elev] = meshgrid(FOV_1,FOV_2);
                FOV_rho = FOV_3*ones(size(FOV_azi));
                [FOV_x,FOV_y,FOV_z] = sph2cart(FOV_azi,FOV_elev,FOV_rho);
                shifted_x = FOV_x - magnet_loc(1);
                shifted_y = FOV_y - magnet_loc(2);
                shifted_z = FOV_z - magnet_loc(3);
                FOV_shifted = ones(numel(FOV_x),3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            end

            angle_rot = -angle;
            FOV_rotated = self.Rotate3D(FOV_shifted,angle_rot,1);

            % Source points of the magnet
            px = [-magnet_dim(1)/2,magnet_dim(1)/2];
            py = [-magnet_dim(2)/2,magnet_dim(2)/2];
            pz = [-magnet_dim(3)/2,magnet_dim(3)/2];
            C = {[1,2],[1,2],[1,2]};
            D = C;
            [D{:}] = ndgrid(C{:});
            zip_list = cell2mat(cellfun(@(m)m(:),D,'uni',0));
            num_elem = length(zip_list);
            % Calculated the field map
            for elem_idx=1:num_elem
                i = zip_list(elem_idx,1);
                j = zip_list(elem_idx,2);
                k = zip_list(elem_idx,3);
                % Distance between the observation point and the source point
                xd = FOV_rotated(:,1)-px(i);
                yd = FOV_rotated(:,2)-py(j);
                zd = FOV_rotated(:,3)-pz(k);
                r = sqrt((xd.^2)+(yd.^2)+(zd.^2));
                
                % Bx
                % For points aligned with the source points, the value
                % inside log() will be zero, but there is always a pair
                % of such value with opposite sign, so we set it to 0
                tmp_sum = zd+r;
                if sum(sum(tmp_sum==0))~=0
                    tmp_sum(tmp_sum==0) = 1;
                end
                Bx_curr = ((-1)^(i+j+k))*log(tmp_sum);
                B_all_local(:,1) = B_all_local(:,1) + Bx_curr;
                % By
                nume = xd.*zd;
                denom = yd.*r;
                B_all_local(:,2) = B_all_local(:,2) + ((-1)^(i+j+k+1))*(atan2(nume,denom));
                % Bz
                tmp_sum = xd+r;
                if sum(sum(tmp_sum==0))~=0
                    tmp_sum(tmp_sum==0) = 1;
                end
                Bz_curr = ((-1)^(i+j+k))*log(tmp_sum);
                B_all_local(:,3) = B_all_local(:,3) + Bz_curr;

            end
            % Multiply Bx/By/Bz by the common factor -Br/4pi
            B_all_local = -1*B_all_local*Br/(4*pi);
            % Change back the coordinates
            angle_rot = angle;
            B_all = self.Rotate3D(B_all_local,angle_rot,2);
            Bx = reshape(B_all(:,1),FOV2_len,FOV1_len);
            By = reshape(B_all(:,2),FOV2_len,FOV1_len);
            Bz = reshape(B_all(:,3),FOV2_len,FOV1_len);
        end

        % Compute the 2D magnetic field produced by multiple cuboid magnets
        function [Bx,By,Bz] = Field2DMultiCuboid(self,FOV1_len,FOV2_len,n_mag,FOV_1,FOV_2,FOV_3,magnet_dim,magnet_loc,angle,Br,surface)
        % INPUT:
        %       FOV1/2_len - the length of FOV_1/2 vector
        %       n_mag - number of magnets
	    %       FOV_1/2[1,m1]/[1,m2] - FOV 1D arrays of two directions
	    %       FOV_3 - a scalar, the third dimension of field calculation
	    %       magnet_dim - the size of all magnets
	    %       magnet_loc - the center locations of all magnets
	    %       angle - the orientation angle of all magnets in the form of (yaw,pitch,row)
	    %       Br - remanence of all magnets
	    %       surface - z: XY plane, y: XZ plane, x: YZ plane, where ab-plane means FOV_1 = FOV_a, FOV_2 = FOV_b
	    % OUTPUT:
	    %       Bx/By/Bz - the x/y/z field components in (T)

            N_FOV = FOV1_len*FOV2_len;
            Bx_local = zeros(n_mag,N_FOV);
            By_local = zeros(n_mag,N_FOV);
            Bz_local = zeros(n_mag,N_FOV);

            % Change the coordinates
            % normal vector = y: XZ slice, z: XY slice, x: YZ slice
            pos_rep_x = reshape(repmat(reshape(magnet_loc(:,1),1,[]),FOV2_len,1),[],1);
            pos_rep_y = reshape(repmat(reshape(magnet_loc(:,2),1,[]),FOV2_len,1),[],1);
            pos_rep_z = reshape(repmat(reshape(magnet_loc(:,3),1,[]),FOV2_len,1),[],1);
            if surface == 'x'
                [FOV_y,FOV_z] = meshgrid(FOV_1,FOV_2);
                shifted_x = ones(FOV2_len*n_mag,FOV1_len)*FOV_3 - pos_rep_x;
                shifted_y = repmat(FOV_y,n_mag,1) - pos_rep_y;
                shifted_z = repmat(FOV_z,n_mag,1) - pos_rep_z;
                FOV_shifted = ones(numel(FOV_y)*n_mag,3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            elseif surface == 'y'
                [FOV_x,FOV_z] = meshgrid(FOV_1,FOV_2);
                shifted_x = repmat(FOV_x,n_mag,1) - pos_rep_x;
                shifted_y = ones(FOV2_len*n_mag,FOV1_len)*FOV_3 - pos_rep_y;
                shifted_z = repmat(FOV_z,n_mag,1) - pos_rep_z;
                FOV_shifted = ones(numel(FOV_x)*n_mag,3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            elseif surface == 'z'
                [FOV_x,FOV_y] = meshgrid(FOV_1,FOV_2);
                shifted_x = repmat(FOV_x,n_mag,1) - pos_rep_x;
                shifted_y = repmat(FOV_y,n_mag,1) - pos_rep_y;
                shifted_z = ones(FOV2_len*n_mag,FOV1_len)*FOV_3 - pos_rep_z;
                FOV_shifted = ones(numel(FOV_x)*n_mag,3);
                FOV_shifted(:,1) = reshape(shifted_x',[],1);
                FOV_shifted(:,2) = reshape(shifted_y',[],1);
                FOV_shifted(:,3) = reshape(shifted_z',[],1);
            else
                % Curvy surface
                [FOV_phi,FOV_z] = meshgrid(FOV_1,FOV_2);
                FOV_rho = FOV_3*ones(size(FOV_phi));
                [FOV_x,FOV_y] = pol2cart(FOV_phi+pi/2,FOV_rho);
                shifted_x = FOV_x - magnet_loc(1);
                shifted_y = FOV_y - magnet_loc(2);
                shifted_z = FOV_z - magnet_loc(3);
                FOV_shifted = ones(numel(FOV_x),3);
                FOV_shifted(:,1) = reshape(shifted_x,[],1);
                FOV_shifted(:,2) = reshape(shifted_y,[],1);
                FOV_shifted(:,3) = reshape(shifted_z,[],1);
            end

            % Rotate the FoV for every magnet
            FOV_rotated = zeros(size(FOV_shifted));
            FOV_final_x = zeros(n_mag,N_FOV);
            FOV_final_y = zeros(n_mag,N_FOV);
            FOV_final_z = zeros(n_mag,N_FOV);
            for mag_idx=1:n_mag
                angle_rot = -angle(mag_idx,:);
                idx_start = (mag_idx-1)*N_FOV + 1;
                idx_end = mag_idx*N_FOV;
                FOV_curr = FOV_shifted(idx_start:idx_end,:);
                FOV_rotated(idx_start:idx_end,:) = self.Rotate3D(FOV_curr,angle_rot,1);
                % Reshape the FoV matrix to separate x/y/z components with shape [n_mag,N_FOV]
                FOV_final_x(mag_idx,:) = reshape(FOV_rotated(idx_start:idx_end,1),1,N_FOV);
			    FOV_final_y(mag_idx,:) = reshape(FOV_rotated(idx_start:idx_end,2),1,N_FOV);
			    FOV_final_z(mag_idx,:) = reshape(FOV_rotated(idx_start:idx_end,3),1,N_FOV);
            end

            % Source points of the magnet
            px = [-magnet_dim(:,1)/2,magnet_dim(:,1)/2];
            py = [-magnet_dim(:,2)/2,magnet_dim(:,2)/2];
            pz = [-magnet_dim(:,3)/2,magnet_dim(:,3)/2];
            p_idx_list = {[1,2],[1,2],[1,2]};
            idx_list_copy = p_idx_list;
            [idx_list_copy{:}] = ndgrid(p_idx_list{:});
            zip_list = cell2mat(cellfun(@(m)m(:),idx_list_copy,'uni',0));
            num_elem = length(zip_list);
            % Calculated the field map
            for elem_idx=1:num_elem
                i = zip_list(elem_idx,1);
                j = zip_list(elem_idx,2);
                k = zip_list(elem_idx,3);
                % Distance between the observation point and the source point
                xd = FOV_final_x-px(:,i);
                yd = FOV_final_y-py(:,j);
                zd = FOV_final_z-pz(:,k);
                r = sqrt((xd.^2)+(yd.^2)+(zd.^2));
                
                % Bx
                % For points aligned with the source points, the value
                % inside log() will be zero, but there is always a pair
                % of such value with opposite sign, so we set it to 0
                tmp_sum = zd+r;
                if sum(sum(tmp_sum==0))~=0
                    tmp_sum(tmp_sum==0) = 1;
                end
                Bx_curr = ((-1)^(i+j+k))*log(tmp_sum);
                Bx_local = Bx_local + Bx_curr;
                % By
                nume = xd.*zd;
                denom = yd.*r;
                By_local = By_local + ((-1)^(i+j+k+1))*(atan2(nume,denom));
                % Bz
                tmp_sum = xd+r;
                if sum(sum(tmp_sum==0))~=0
                    tmp_sum(tmp_sum==0) = 1;
                end
                Bz_curr = ((-1)^(i+j+k))*log(tmp_sum);
                Bz_local = Bz_local + Bz_curr;

            end
            % Change back the coordinates for every magnet
            Bx_all = zeros(size(Bx_local));
            By_all = zeros(size(Bx_local));
            Bz_all = zeros(size(Bx_local));
            for mag_idx=1:n_mag
                angle_rot = angle(mag_idx,:);
                B_all_curr = [Bx_local(mag_idx,:);By_local(mag_idx,:);Bz_local(mag_idx,:)];
                B_all_rot = self.Rotate3D(B_all_curr',angle_rot,2);
                Bx_all(mag_idx,:) = B_all_rot(:,1);
			    By_all(mag_idx,:) = B_all_rot(:,2);
			    Bz_all(mag_idx,:) = B_all_rot(:,3);
            end
            
            % Multiply Bx/By/Bz by the common factor -Br/4pi
            Bx_all = -1*Br*Bx_all/(4*pi);
            By_all = -1*Br*By_all/(4*pi);
            Bz_all = -1*Br*Bz_all/(4*pi);
            Bx = reshape(Bx_all,FOV1_len,FOV2_len)';
            By = reshape(By_all,FOV1_len,FOV2_len)';
            Bz = reshape(Bz_all,FOV1_len,FOV2_len)';
        end
        
        % ========= Magnetic Force Calculation ========== %
        % Static force acting on a single magnet from arbitrary number of sources
        function F_target = ForceSingle(self,target_idx,source_idx,FOV_del_ratio)
        % INPUT:
        %       target_idx - In the form of [type,idx],
        %                    where type: 1 - cuboid, 2 - fan
        %       source_idx[n,2] - Each is in the form of [type,idx]
        %       FOV_del_ratio - The resolution of division for each side
        % OUTPUT:
        %       F_target[1,3] - The x, y & z components of force acting on the
        %                       target magnet
           
            
            % FOV creation
            % 1 - cuboid magnet: x-z surfaces
            if target_idx(1) == 1
                surface = 'y';
                mag_dim_target = self.cub_dim(target_idx(2),:);
                FOV_1 = -mag_dim_target(1)/2:(mag_dim_target(1)/FOV_del_ratio):mag_dim_target(1)/2;
                FOV_2 = -mag_dim_target(3)/2:(mag_dim_target(3)/FOV_del_ratio):mag_dim_target(3)/2;
                FOV_3_list = [-mag_dim_target(2)/2,mag_dim_target(2)/2];
                FOV_del_S = mag_dim_target(1)*mag_dim_target(3)/FOV_del_ratio.^2;
                % Load other parameters
                loc_target = self.cub_loc(target_idx(2),:);
                angle_target = self.cub_angle(target_idx(2),:);
                Br_target = self.cub_Br(target_idx(2));
                pol_target = 1;
            end
            
            % ===== Coordinate transformation =====
            % Load source locations (cuboid: center loc, fan: local origin)
            loc_other_raw = zeros(size(source_idx,1),3);
            angle_other_raw = zeros(size(source_idx,1),3);
            for magnet_idx=1:size(source_idx,1)
                if source_idx(magnet_idx,1) == 1
                    loc_other_raw(magnet_idx,:) = self.cub_loc(source_idx(magnet_idx,2),:);
                    angle_other_raw(magnet_idx,:) = self.cub_angle(source_idx(magnet_idx,2),:);
                end
            end
            % Move to the new origin
            loc_other_unrot = loc_other_raw - loc_target;
            angle_other = angle_other_raw - angle_target;
            % Rotate the coordinates by angle_target
            angle_rot = angle_target;
            loc_other = self.Rotate3D(loc_other_unrot,-angle_rot,1);

            % Force calculation
            F_target = zeros(1,3);
            mu_air = 1.25663753e-6;
            sigma_m = pol_target*Br_target/mu_air;
            
            % Loop through each magnet
            for magnet_idx=1:size(source_idx,1)
                FOV_3_counter = 1;
                for FOV_3=FOV_3_list
                    if target_idx(1) == 1
                        % Surface normal vector direction
                        surface_n = sign(FOV_3);
                    else
                        % Surface normal vector direction
                        surface_n = (-1)^FOV_3_counter;
                        % For fan magnet as target, get the current surface FOV
                        if FOV_3_counter == 1
                            %FOV_1 = FOV_1_in;
                            FOV_del_S = FOV_del_S_1;
                        else
                            %FOV_1 = FOV_1_out;
                            FOV_del_S = FOV_del_S_2;
                        end
                    end

                    loc_source = loc_other(magnet_idx,:);
                    angle_source = angle_other(magnet_idx,:);
                    % Cuboid or fan-shaped source magnet
                    if source_idx(magnet_idx,1) == 1        % Cuboid
                        % Calculate the magnetic field on the surfaces
                        dim_source = self.cub_dim(source_idx(magnet_idx,2),:);
                        Br_source = self.cub_Br(source_idx(magnet_idx,2));
                        [Bx,By,Bz] = self.Field2DSingleCuboid(length(FOV_1),length(FOV_2),FOV_1,FOV_2,FOV_3,dim_source,loc_source,angle_source,Br_source,surface);
                    end
                    
                    F_target(1) = F_target(1) + surface_n*sum(Bx,'all')*sigma_m*FOV_del_S*1e-6;
                    F_target(2) = F_target(2) + surface_n*sum(By,'all')*sigma_m*FOV_del_S*1e-6;
                    F_target(3) = F_target(3) + surface_n*sum(Bz,'all')*sigma_m*FOV_del_S*1e-6;

                    FOV_3_counter = FOV_3_counter + 1;
                end
            end
            % Coordinate transformation
            angle_rot = angle_target;
            F_target = self.Rotate3D(F_target,angle_rot,2);
        end
    
        % Static force acting on all magnets from all the other magnets
        function F_all_static = ForceAllStatic(self,FOV_del_ratio)
        % INPUT:
        %       FOV_del_ratio - The resolution of division for each side
        % OUTPUT:
        %       F_target[n,3] - The x, y & z components of force acting on
        %                       all target magnets

            % Loop through all magnets
            cub_idx_list = [ones(self.cub_numMag,1),(1:self.cub_numMag)'];

            total_idx_list = cub_idx_list;
            total_numMag = self.cub_numMag;
            F_all_static = zeros(total_numMag,3);
            for mag_idx=1:total_numMag
                curr_mag_idx = total_idx_list(mag_idx,:);
                other_mag_idx = total_idx_list(setdiff(1:total_numMag,mag_idx),:);
                F_all_static(mag_idx,:) = self.ForceSingle(curr_mag_idx,other_mag_idx,FOV_del_ratio);
                if mod(mag_idx,10) == 0 || mag_idx == total_numMag
                    fprintf('Magnet #%d\n',mag_idx);
                end
            end
        end
    
        % Static force acting on multiple magnets from arbitrary number of sources
        function F_target_list = ForceFree(self,target_idx_list,source_idx_list,FOV_del_ratio)
        % INPUT:
        %       target_idx_list[n,2] - In the form of [type,idx],
        %                    where type: 1 - cuboid, 2 - fan
        %       source_idx[n,2] - Each is in the form of [type,idx]
        %       FOV_del_ratio - The resolution of division for each side
        % OUTPUT:
        %       F_target[1,3] - The x, y & z components of force acting on the
        %                       target magnet

            % Loop through all magnets in the given list
            total_numMag = size(target_idx_list,1);
            F_target_list = zeros(total_numMag,3);
            for mag_idx=1:total_numMag
                curr_mag_idx = target_idx_list(mag_idx,:);
                other_mag_idx = source_idx_list;
                F_target_list(mag_idx,:) = self.ForceSingle(curr_mag_idx,other_mag_idx,FOV_del_ratio);
                if mod(mag_idx,10) == 0 || mag_idx == total_numMag
                    fprintf('Magnet #%d\n',mag_idx);
                end
            end
        end

        % Static force acting on multiple magnets IN A GROUP from arbitrary number of sources
        function F_target = ForceGroup(self,target_idx_list,source_idx_list,FOV_del_ratio)
        % INPUT:
        %       target_idx_list[n,2] - In the form of [type,idx],
        %                    where type: 1 - cuboid, 2 - fan
        %       source_idx[n,2] - Each is in the form of [type,idx]
        %       FOV_del_ratio - The resolution of division for each side
        % OUTPUT:
        %       F_target[1,3] - The x, y & z components of force acting on the
        %                       target magnet

            % Loop through all magnets in the given list
            total_numMag = size(target_idx_list,1);
            F_target_list = zeros(total_numMag,3);
            for mag_idx=1:total_numMag
                curr_mag_idx = target_idx_list(mag_idx,:);
                other_mag_idx = source_idx_list;
                F_target_list(mag_idx,:) = self.ForceSingle(curr_mag_idx,other_mag_idx,FOV_del_ratio);
                if mod(mag_idx,10) == 0 || mag_idx == total_numMag
                    fprintf('Magnet #%d\n',mag_idx);
                end
            end
            F_target = sum(F_target_list,1);
        end
    
        % Estimate the total weight of PMA (magnet only) in [kg]
        function weight_total = Weight(self)
        % OUTPUT:
        %       weight_total - the total weight of the magnets

            % Reference: A 10x10x10mm^3 magnet weights 8g
            weight_ref = 8e-3;         % [kg]
            volume_ref = 1e3;          % [mm^3]
            % Cuboid
            volume_cub = 0;
            if self.cub_numMag > 0
                volume_cub = sum(prod(self.cub_dim,2));
            end
            weight_total = volume_cub/volume_ref*weight_ref;
        end
    end

    methods (Static)
        % Perform 3D rotation on (1) location w.r.t. the origin or (2) 3D force or (3) 3D magnetic field
        function var_final = Rotate3D(var_orig,angle_rot,mode)
        % INPUT:
        %       var_orig - the original quantity to be rotated
        %       angle_rot - the rotation angle in the form of (yaw,pitch,roll) and the unit of [degree]
        %       mode - the rotation order with 1: roll->pitch->yaw, and 2: yaw->pitch->roll
            angle_yaw = angle_rot(1);
            angle_pitch = angle_rot(2);
            angle_roll = angle_rot(3);
            rot_roll = [1,0,0;...
                        0,cosd(angle_roll),-sind(angle_roll);...
                        0,sind(angle_roll),cosd(angle_roll)];
            rot_pitch = [cosd(angle_pitch),0,sind(angle_pitch);...
                         0,1,0;...
                         -sind(angle_pitch),0,cosd(angle_pitch)];
            rot_yaw = [cosd(angle_yaw),-sind(angle_yaw),0;...
                       sind(angle_yaw),cosd(angle_yaw),0;...
                       0,0,1];
            if mode == 1
                rot_mat = (rot_roll*rot_pitch*rot_yaw)';
            elseif mode == 2
                rot_mat = (rot_yaw*rot_pitch*rot_roll)';
            end
            var_final = var_orig*rot_mat;
        end
    end
end