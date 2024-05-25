function transducer_object = PFC_Make3DTransducer(transducer_info)
% Make 3D mesh transducer using 'delaunayTriangulation'
% Using the 'right-handed' Cartesian coordinate system
% The basic normal vector should be (0, 0, 1) direction
% Subfunctions : this_InnerVerteciesHex2D, this_RotationMatrix
%
% Basic Inputs : Frequency[Hz], Pressure[Pa], Phase[rad], Type
% Additional Inputs : rotate[degree], position[m], medium, max_edgelength[m]
%
% transducer_info.Frequency = 1e6;
% transducer_info.Pressure = 1;
% transducer_info.Phase = 0;
%
% transducer_info.rotate = [0, 0, 0];
% transducer_info.position = [0, 0, 0];
% transducer_info.max_edgelength = 0.2e-3;
%
%
% Transdycer Type
% Basic Type : 'circle', 'ring', 'ring_angular', 'square', 'rect'
%
% transducer_info.Type = 'circle';
% transducer_info.Radius = 5e-3;
%
% transducer_info.Type = 'ring';
% transducer_info.Radius_in = 5e-3;
% transducer_info.Radius_out = 6e-3;
%
% transducer_info.Type = 'ring_angular';
% transducer_info.Radius_in = 5e-3;
% transducer_info.Radius_out = 6e-3;
% transducer_info.Ring_angle = 315; % deg
%
% transducer_info.Type = 'square';
% transducer_info.Length = 5e-3;
%
% transducer_info.Type = 'rect';
% transducer_info.Length_x = 5e-3;
% transducer_info.Length_y = 10e-3;
%
%
% Array Type : 'PMUT_MUTi_2.1', 'Sparse_CMUT_v1_cell', 'Sparse_CMUT_v1'
%
% transducer_info.Type = 'PMUT_MUTi_2.1';
% transducer_info.focal_length = 12e-3;
%
% transducer_info.Type = 'Sparse_CMUT_v1_cell';
%
% transducer_info.Type = 'Sparse_CMUT_v1';
% transducer_info.chan_list = 1:128;
%
% transducer_object = PFC_Make3DTransducer(transducer_info);
%
%
% Projection Type : 'spherical_cap'
%
% transducer_info.Type = 'spherical_cap';
% transducer_info.Radius = 10e-3;
% transducer_info.Focal_length = 30e-3;
%
% transducer_object = PFC_Make3DTransducer(transducer_info);


%% Basic Input

% Essential Input
Frequency_designed = transducer_info.Frequency;     % [Hz]
Pressure_input = transducer_info.Pressure;          % [Pa]
Phase_input = transducer_info.Phase;                % [rad]

% Additional Input
if(isfield(transducer_info, 'rotate'))
    Transducer_drotate = transducer_info.rotate;        % [degree]
else
    Transducer_drotate = [0, 0, 0];
end


if(isfield(transducer_info, 'position'))
    Transducer_position = transducer_info.position;     % [m]
else
    Transducer_position = [0, 0, 0];
end


if(isfield(transducer_info, 'medium'))
    Medium_name = transducer_info.medium;
else
    Medium_name = 'water';
end

Medium_data = PFC_GetMedium(Medium_name);

Medium_data.wavelength = Medium_data.c/Frequency_designed; % [m]


if(isfield(transducer_info, 'max_edgelength'))
    max_edgelength = transducer_info.max_edgelength; % [m]
else
    max_edgelength = Medium_data.wavelength/10;
end


% Output
transducer_object.FrequencyDesigned = Frequency_designed;
transducer_object.MediumDesigned = Medium_data;


%% transducer Type
% Each type should calculate
% Tri_Points, Tri_Constraints, Tri_ConnectivityList,
% transducer_pressure_gain, transducer_phase_delay

switch(transducer_info.Type)
    case {'circle', 'ring', 'ring_angular', 'square', 'rect'}
        [Tri_Points, Tri_Constraints, Tri_ConnectivityList]...
            = this_Make2DTransducerBasic(transducer_info, max_edgelength);

        transducer_pressure_gain = ones(size(Tri_ConnectivityList, 1), 1);
        transducer_phase_delay = zeros(size(Tri_ConnectivityList, 1), 1);


    case {'spherical_cap'}
        [Tri_Points, Tri_Constraints, Tri_ConnectivityList,...
            transducer_pressure_gain, transducer_phase_delay]...
            = this_Make3DTransducer_Projection(transducer_info, max_edgelength);


    case {'PMUT_MUTi_2.1', 'Sparse_CMUT_v1_cell', 'Sparse_CMUT_v1', 'CMUT_HG'}
        [Tri_Points, Tri_Constraints, Tri_ConnectivityList,...
            transducer_pressure_gain, transducer_phase_delay]...
            = this_Make2DTransducerArray(transducer_info, max_edgelength, Medium_data);
end


%% Rotate
Rxyz = this_RotationMatrix(Transducer_drotate);

Rxyz = Rxyz';

Tri_Points = Tri_Points*Rxyz;


%% Move Position
Tri_Points = Tri_Points + ones(size(Tri_Points, 1), 1)*Transducer_position;


%% Calculate Pressure
transducer_pressure =...
    Pressure_input .* transducer_pressure_gain...
    .* exp(-1i.*Phase_input) .* exp(-1i.*transducer_phase_delay);

transducer_object.Pressure = transducer_pressure;

transducer_object.Tri_Points = Tri_Points;
transducer_object.Tri_ConnectivityList = Tri_ConnectivityList;
transducer_object.Tri_Constraints = Tri_Constraints;

end









%% Other Functions

function [Tri_Points, Tri_Constraints, Tri_ConnectivityList]...
    = this_Make2DTransducerBasic(transducer_info, max_edgelength)

constraints_edge = [];

switch(transducer_info.Type)
    case 'circle'
        circle_radius = transducer_info.Radius;

        % Find Edge Vertex
        dx = max_edgelength/sqrt(2);

        if(isfield(transducer_info, 'dtheta'))
            dtheta = transducer_info.dtheta;
        else
            dtheta = 2*asin(dx/2/circle_radius);
        end

        dtheta_i = ceil(2*pi/dtheta);
        dtheta = 2*pi/dtheta_i;

        theta_list = (1:dtheta_i).*dtheta;
        theta_list = theta_list';
        edge_vertex = circle_radius.*[cos(theta_list), sin(theta_list)];

        % Find Inner Vertex
        inner_vertex = this_InnerVerteciesHex2D(circle_radius, dx);

        inner_vertex_r = sqrt(sum(inner_vertex.^2, 2));

        inner_vertex_on = inner_vertex_r < (circle_radius - 0.2.*dx);

        inner_vertex = inner_vertex(inner_vertex_on, :);

    case 'ring'
        ring_radius_in = transducer_info.Radius_in;
        ring_radius_out = transducer_info.Radius_out;


        % Find Edge Vertex
        dx = max_edgelength/sqrt(2);
        dtheta = 2*asin(dx/2/ring_radius_in);
        dtheta_i = ceil(2*pi/dtheta);
        dtheta = 2*pi/dtheta_i;

        theta_list = (1:dtheta_i).*dtheta;
        theta_list = theta_list';
        inner_edge_vertex = ring_radius_in.*[cos(theta_list), sin(theta_list)];

        dtheta = 2*asin(dx/2/ring_radius_out);
        dtheta_i = ceil(2*pi/dtheta);
        dtheta = 2*pi/dtheta_i;

        theta_list = (1:dtheta_i).*dtheta;
        theta_list = theta_list';
        outer_edge_vertex = ring_radius_out.*[cos(theta_list), sin(theta_list)];

        edge_vertex = [inner_edge_vertex; outer_edge_vertex];


        constraints_in =...
            [(1:(size(inner_edge_vertex, 1) - 1))', (1:(size(inner_edge_vertex, 1) - 1))' + 1;...
            size(inner_edge_vertex, 1), 1];

        constraints_out = size(inner_edge_vertex, 1) +...
            [(1:(size(outer_edge_vertex, 1) - 1))', (1:(size(outer_edge_vertex, 1) - 1))' + 1;...
            size(outer_edge_vertex, 1), 1];

        constraints_edge = [constraints_in; constraints_out];


        % Find Inner Vertex
        inner_vertex = this_InnerVerteciesHex2D(ring_radius_out, dx);

        inner_vertex_r = sqrt(sum(inner_vertex.^2, 2));

        inner_vertex_on = (inner_vertex_r > (ring_radius_in + 0.2.*dx)) &...
            (inner_vertex_r < (ring_radius_out - 0.2.*dx));

        inner_vertex = inner_vertex(inner_vertex_on, :);

    case 'ring_angular'
        ring_radius_in = transducer_info.Radius_in;
        ring_radius_out = transducer_info.Radius_out;
        ring_angle = transducer_info.Ring_angle;
        ring_angle_rad = pi/180*ring_angle;


        % Find Edge Vertex
        dx = max_edgelength/sqrt(2);
        dtheta = 2*asin(dx/2/ring_radius_in);
        dtheta_i = ceil(ring_angle_rad/dtheta);
        dtheta = ring_angle_rad/dtheta_i;

        theta_list = (0:dtheta_i).*dtheta;
        theta_list = theta_list';
        inner_edge_vertex = ring_radius_in.*[cos(theta_list), sin(theta_list)];


        dx = max_edgelength/sqrt(2);
        dtheta = 2*asin(dx/2/ring_radius_out);
        dtheta_i = ceil(ring_angle_rad/dtheta);
        dtheta = ring_angle_rad/dtheta_i;

        theta_list = (0:dtheta_i).*dtheta;
        theta_list = theta_list';
        outer_edge_vertex = ring_radius_out.*[cos(theta_list), sin(theta_list)];

        edge_vertex_radius = (ring_radius_in + dx/2):dx:(ring_radius_out - dx/2);

        zero_deg_edge_vertex = edge_vertex_radius'.*[cos(0), sin(0)];
        angle_deg_edge_vertex = edge_vertex_radius'.*[cos(ring_angle_rad), sin(ring_angle_rad)];

        edge_vertex = [zero_deg_edge_vertex; outer_edge_vertex;...
            flipud(angle_deg_edge_vertex); flipud(inner_edge_vertex)];

        constraints_edge =...
            [(1:(size(edge_vertex, 1) - 1))', (1:(size(edge_vertex, 1) - 1))' + 1;...
            size(edge_vertex, 1), 1];


        % Find Inner Vertex
        inner_vertex = this_InnerVerteciesHex2D(ring_radius_out, dx);

        inner_vertex_r = sqrt(sum(inner_vertex.^2, 2));

        inner_vertex_a = mod(atan2(inner_vertex(:, 2), inner_vertex(:, 1)), 2*pi);

        inner_vertex_on =...
            (inner_vertex_r > (ring_radius_in + 0.2.*dx)) &...
            (inner_vertex_r < (ring_radius_out - 0.2.*dx)) &...
            (inner_vertex_a < (ring_angle_rad - dtheta)) &...
            (inner_vertex_a > dtheta);

        inner_vertex = inner_vertex(inner_vertex_on, :);

    case 'square'
        square_length = transducer_info.Length;
        square_length_half = square_length/2;


        % Find Edge Vertex
        dx = max_edgelength/sqrt(2);
        dx_i = ceil(square_length_half/dx);
        dx = square_length_half/dx_i;

        x_list = ((-dx_i + 1):dx_i).*dx;
        x_list = x_list';
        y_list = square_length_half.*ones(size(x_list));

        edge_vertex = [...
            x_list, y_list;...
            -x_list, -y_list;...
            y_list, -x_list;...
            -y_list, x_list];


        % Find Inner Vertex
        inner_vertex = this_InnerVerteciesHex2D(square_length, dx);

        inner_vertex_on = (abs(inner_vertex(:, 1)) < (square_length_half - 0.2.*dx)) &...
            (abs(inner_vertex(:, 2)) < (square_length_half - 0.2.*dx));

        inner_vertex = inner_vertex(inner_vertex_on, :);

    case 'rect'
        rect_length_x = transducer_info.Length_x;
        rect_length_y = transducer_info.Length_y;

        rect_length_x_half = rect_length_x/2;
        rect_length_y_half = rect_length_y/2;


        % Find Edge Vertex
        dxx = max_edgelength/sqrt(2);
        dxx_i = ceil(rect_length_x_half/dxx);
        dxx = rect_length_x_half/dxx_i;

        xx_list = ((-dxx_i + 1):dxx_i).*dxx;
        xx_list = xx_list';
        xy_list = rect_length_y_half.*ones(size(xx_list));

        dxy = max_edgelength/sqrt(2);
        dxy_i = ceil(rect_length_y_half/dxy);
        dxy = rect_length_y_half/dxy_i;

        yy_list = ((-dxy_i + 1):dxy_i).*dxy;
        yy_list = yy_list';
        yx_list = rect_length_x_half.*ones(size(yy_list));


        edge_vertex = [...
            xx_list, xy_list;...
            -xx_list, -xy_list;...
            yx_list, -yy_list;...
            -yx_list, yy_list];

        dx = min(dxx, dxy);


        % Find Inner Vertex
        inner_vertex = this_InnerVerteciesHex2D(max(rect_length_x, rect_length_y), dx);

        inner_vertex_on = (abs(inner_vertex(:, 1)) < (rect_length_x_half - 0.2.*dx)) &...
            (abs(inner_vertex(:, 2)) < (rect_length_y_half - 0.2.*dx));

        inner_vertex = inner_vertex(inner_vertex_on, :);
end


% Add All Vertex
total_vertex = [edge_vertex; inner_vertex];

% Delaunay Triangulation

if(isempty(constraints_edge))
    DT = delaunayTriangulation(total_vertex);
    Tri_ConnectivityList = DT.ConnectivityList;
else
    DT = delaunayTriangulation(total_vertex, constraints_edge);
    TF = isInterior(DT);
    Tri_ConnectivityList = DT.ConnectivityList(TF, :);
end

% Calculated Data
Tri_Points = [DT.Points, zeros(size(DT.Points, 1), 1)];
Tri_Constraints = DT.Constraints;

end


function [Tri_Points, Tri_Constraints, Tri_ConnectivityList,...
    transducer_pressure_gain, transducer_phase_delay]...
    = this_Make3DTransducer_Projection(transducer_info, max_edgelength)

switch(transducer_info.Type)
    case 'spherical_cap'
        sc_radius = transducer_info.Radius;
        sc_focal_length = transducer_info.Focal_length;

        % Find Edge Vertex
        dx = max_edgelength/sqrt(2);

        dtheta = 2*asin(dx/2/sc_radius);

        dtheta_i = ceil(2*pi/dtheta);
        dtheta = 2*pi/dtheta_i;

        theta_list = (1:dtheta_i).*dtheta;
        theta_list = theta_list';
        edge_vertex = sc_radius.*[cos(theta_list), sin(theta_list)];

        % Find Inner Vertex
        inner_vertex = this_InnerVerteciesHex2D(sc_radius, dx);

        inner_vertex_r = sqrt(sum(inner_vertex.^2, 2));

        inner_vertex_on = inner_vertex_r < (sc_radius - 0.2.*dx);

        inner_vertex = inner_vertex(inner_vertex_on, :);

        % Add All Vertex
        total_vertex_proj = [edge_vertex; inner_vertex];

        % Delaunay Triangulation
        DT = delaunayTriangulation(total_vertex_proj);
        Tri_ConnectivityList = DT.ConnectivityList;

        % Calculated Data
        Tri_Point_z = sc_focal_length - sqrt(sc_focal_length.^2 - sum(DT.Points.^2, 2));

        Tri_Points = [DT.Points, Tri_Point_z];
        Tri_Constraints = DT.Constraints;

        transducer_pressure_gain = ones(size(Tri_ConnectivityList, 1), 1);
        transducer_phase_delay = zeros(size(Tri_ConnectivityList, 1), 1);
end

end


function [Tri_Points, Tri_Constraints, Tri_ConnectivityList,...
    transducer_pressure_gain, transducer_phase_delay]...
    = this_Make2DTransducerArray(transducer_info, max_edgelength, medium_data)

switch(transducer_info.Type)
    case 'PMUT_MUTi_2.1'
        fl = transducer_info.focal_length;

        r_in = [...
            1.04e-3,...
            3.28e-3,...
            4.42e-3,...
            5.56e-3,...
            6.70e-3,...
            7.58e-3,...
            8.46e-3,...
            9.34e-3]./2;

        r_out = [...
            2.6e-3,...
            4.09e-3,...
            5.23e-3,...
            6.37e-3,...
            7.25e-3,...
            8.13e-3,...
            9.01e-3,...
            9.89e-3]./2;

        r_center = (r_in + r_out)/2;

        fl_diff = fl - sqrt(fl.^2 + r_center.^2);

        Phase_delay_list = 2*pi*fl_diff/medium_data.wavelength;

        transducer_info_i = transducer_info;
        transducer_info_i.Pressure = 1;

        transducer_info_i.rotate = [0, 0, 0];
        transducer_info_i.position = [0, 0, 0];

        transducer_info_i.Type = 'ring_angular';
        transducer_info_i.max_edgelength = max_edgelength;

        transducer_array = [];

        for i = 1:length(r_center)
            transducer_info_i.Radius_in = r_in(i);
            transducer_info_i.Radius_out = r_out(i);
            transducer_info_i.Ring_angle = 315; % deg

            transducer_info_i.Phase = Phase_delay_list(i);

            transducer_object_i = PFC_Make3DTransducer(transducer_info_i);

            transducer_array = this_add(transducer_array, transducer_object_i);
        end

        transducer_object_i = PFC_TransducerArray2Transducer(transducer_array);

        transducer_pressure = transducer_object_i.Pressure;

        Tri_Points = transducer_object_i.Tri_Points;
        Tri_Constraints = transducer_object_i.Tri_Constraints;
        Tri_ConnectivityList = transducer_object_i.Tri_ConnectivityList;

        transducer_pressure_gain = transducer_pressure;
        transducer_phase_delay = zeros(size(transducer_pressure));

    case 'Sparse_CMUT_v1_cell'
        CMUT_cell_r = 24e-6;
        CMUT_cell_d = 53e-6;

        transducer_info_i = transducer_info;
        transducer_info_i.Pressure = 1;

        transducer_info_i.rotate = [0, 0, 0];
        transducer_info_i.position = [0, 0, 0];

        transducer_info_i.Type = 'circle';
        transducer_info_i.max_edgelength = max_edgelength;

        transducer_info_i.Radius = CMUT_cell_r;

        transducer_array = [];

        for i = -1:1
            for j = -1:1
                transducer_info_i.position = [CMUT_cell_d.*i, CMUT_cell_d.*j, 0];

                transducer_object_i = PFC_Make3DTransducer(transducer_info_i);

                transducer_array = this_add(transducer_array, transducer_object_i);
            end
        end

        transducer_object_i = PFC_TransducerArray2Transducer(transducer_array);

        transducer_pressure = transducer_object_i.Pressure;

        Tri_Points = transducer_object_i.Tri_Points;
        Tri_Constraints = transducer_object_i.Tri_Constraints;
        Tri_ConnectivityList = transducer_object_i.Tri_ConnectivityList;

        transducer_pressure_gain = transducer_pressure;
        transducer_phase_delay = zeros(size(transducer_pressure));

    case 'Sparse_CMUT_v1'
        chan_list = transducer_info.chan_list;

        device_pos_001 = [...
            12.5, -11.5;...
            15.5, -14.5;...
            14.5, -13.5;...
            13.5, -10.5;...
            12.5, -8.5;...
            15.5, -10.5;...
            7.5, -5.5;...
            9.5, -5.5;...
            11.5, -5.5;...
            10.5, -4.5;...
            13.5, -5.5;...
            13.5, -4.5;...
            9.5, -2.5;...
            10.5, -2.5;...
            6.5, -2.5;...
            15.5, -0.5];

        device_pos_002 = flipud(device_pos_001) .* [ones(16,1), -ones(16,1)];

        device_pos_003 = [...
            11.5, 10.5;...
            6.5, 5.5;...
            13.5, 15.5;...
            12.5, 13.5;...
            11.5, 13.5;...
            7.5, 7.5;...
            7.5, 10.5;...
            8.5, 12.5;...
            6.5, 12.5;...
            6.5, 13.5;...
            6.5, 10.5;...
            2.5, 3.5;...
            2.5, 4.5;...
            2.5, 6.5;...
            1.5, 15.5;...
            0.5, 6.5];

        device_pos_004 = flipud(device_pos_003) .* [-ones(16,1), ones(16,1)];

        device_pos_1 = [device_pos_001; device_pos_002; device_pos_003; device_pos_004];
        device_pos_2 = device_pos_1 .* [-ones(64,1), -ones(64,1)];

        device_pos = [device_pos_1; device_pos_2];

        device_pos = device_pos * 180e-6;

        device_pos(:,3) = 0;


        transducer_info_i = transducer_info;
        transducer_info_i.Pressure = 1;

        transducer_info_i.rotate = [0, 0, 0];
        transducer_info_i.position = [0, 0, 0];

        transducer_info_i.Type = 'Sparse_CMUT_v1_cell';
        transducer_info_i.max_edgelength = max_edgelength;

        transducer_array = [];

        for i = 1:length(chan_list)
            chan_i = chan_list(i);

            transducer_info_i.position = device_pos(chan_i, :);

            transducer_object_i = PFC_Make3DTransducer(transducer_info_i);

            transducer_array = this_add(transducer_array, transducer_object_i);
        end

        transducer_object_i = PFC_TransducerArray2Transducer(transducer_array);

        transducer_pressure = transducer_object_i.Pressure;

        Tri_Points = transducer_object_i.Tri_Points;
        Tri_Constraints = transducer_object_i.Tri_Constraints;
        Tri_ConnectivityList = transducer_object_i.Tri_ConnectivityList;

        transducer_pressure_gain = transducer_pressure;
        transducer_phase_delay = zeros(size(transducer_pressure));
    case 'CMUT_HG'

                 fl = transducer_info.focal_length;

        r_in = [...
            0.015e-3,...
            0.7805e-3,...
            1.3513e-3,...
            1.7472e-3,...
            2.0711e-3,...
            2.3532e-3,...
            2.607e-3,...
            2.8401e-3,...
            3.0573e-3,...
            3.2123e-3,...
            3.3673e-3];

        r_out = [...
            0.7705e-3,...
            1.3413e-3,...
            1.7372e-3,...
            2.0611e-3,...
            2.3432e-3,...
            2.597e-3,...
            2.8301e-3,...
            3.0473e-3,...
            3.2023e-3,...
            3.3573e-3,...
            3.5123e-3];


            a = 0; % x Location
            b = 0; % y Location


        r_center = (r_in + r_out)/2;

                
                r_diff = sqrt((r_center - a).^2 + b^2);

                fl_diff = fl - sqrt(fl^2 + r_diff.^2);

                Phase_delay_list = 2*pi*fl_diff/medium_data.wavelength;
        % 1 sector
        transducer_info_i = transducer_info;
        transducer_info_i.Pressure = 1;
        transducer_info_i.rotate = [0, 0, 0];
        transducer_info_i.position = [0, 0, 0];

        transducer_info_i.Type = 'ring_angular';
        transducer_info_i.max_edgelength = max_edgelength;

        transducer_array = [];

        for i = 1:length(r_center)
            transducer_info_i.Radius_in = r_in(i);
            transducer_info_i.Radius_out = r_out(i);
            transducer_info_i.Ring_angle = 180; % deg

            transducer_info_i.Phase = Phase_delay_list(i);
%             transducer_info_i.Phase = 0;

            transducer_object_i = PFC_Make3DTransducer(transducer_info_i);

            transducer_array = this_add(transducer_array, transducer_object_i);
        end


                r_diff = sqrt((r_center + a).^2 + b^2);

                fl_diff = fl - sqrt(fl^2 + r_diff.^2);

                Phase_delay_list = 2*pi*fl_diff/medium_data.wavelength + pi;

        % 2 sector
        transducer_info_ii = transducer_info;
        transducer_info_ii.Pressure = 1;
        transducer_info_ii.rotate = [0, 0, 180];
        transducer_info_ii.position = [0, 0, 0];

        transducer_info_ii.Type = 'ring_angular';
        transducer_info_ii.max_edgelength = max_edgelength;


        for i = 1:length(r_center)
            transducer_info_ii.Radius_in = r_in(i);
            transducer_info_ii.Radius_out = r_out(i);
            transducer_info_ii.Ring_angle = 180; % deg

            transducer_info_ii.Phase = Phase_delay_list(i);
%             transducer_info_ii.Phase = 0;

            transducer_object_ii = PFC_Make3DTransducer(transducer_info_ii);

            transducer_array = this_add(transducer_array, transducer_object_ii);
        end


        % transducer_info.Frequency = 5e6;
        % transducer_info.Pressure = 1;
        % transducer_info.Phase = 2*pi*(fl - sqrt(fl^2 + fl_xy^2))/medium_data.wavelength;
        % transducer_info.position = [0, 0, 0];
        % transducer_info.Type = 'circle';
        % transducer_info.Radius = 0.544e-3;
        % 
        % transducer_object_i = PFC_Make3DTransducer(transducer_info);
        % transducer_array = this_add(transducer_array, transducer_object_i);
%         %
        transducer_object_i = PFC_TransducerArray2Transducer(transducer_array);

        transducer_pressure = transducer_object_i.Pressure;

        Tri_Points = transducer_object_i.Tri_Points;
        Tri_Constraints = transducer_object_i.Tri_Constraints;
        Tri_ConnectivityList = transducer_object_i.Tri_ConnectivityList;

        transducer_pressure_gain = transducer_pressure;
        transducer_phase_delay = zeros(size(transducer_pressure));
end

end


function inner_vertex = this_InnerVerteciesHex2D(max_length, dx)

v1 = dx.*[1, 0];
v2 = dx.*[cosd(60), sind(60)];
v3 = dx.*[cosd(120), sind(120)];

dx_grid_i = ceil(max_length/dx);

v_i = 1:dx_grid_i;
v_ii = 0:dx_grid_i;

[v_i_mesh, v_ii_mesh] = meshgrid(v_i, v_ii);

v_x = reshape(v_i_mesh, [numel(v_i_mesh) ,1]);
v_y = reshape(v_ii_mesh, [numel(v_ii_mesh) ,1]);

v12 = v_x*v1 + v_y*v2;
v23 = v_x*v2 + v_y*v3;
v31 = v_x*v3 - v_y*v1;

inner_vertex_12 = [v12; -v12];
inner_vertex_23 = [v23; -v23];
inner_vertex_31 = [v31; -v31];

inner_vertex = [[0, 0]; inner_vertex_12; inner_vertex_23; inner_vertex_31];

end


function Rxyz = this_RotationMatrix(drotate)

x_rot = drotate(1);
y_rot = drotate(2);
z_rot = drotate(3);

Rx = [...
    1, 0, 0;...
    0, cosd(x_rot), -sind(x_rot);
    0, sind(x_rot), cosd(x_rot)];

Ry = [...
    cosd(y_rot), 0, sind(y_rot);...
    0, 1, 0;
    -sind(y_rot), 0, cosd(y_rot)];

Rz = [...
    cosd(z_rot), -sind(z_rot), 0;
    sind(z_rot), cosd(z_rot), 0;
    0, 0, 1];

Rxyz = Rz*Ry*Rx;

end


function list_o = this_add(list_i, add_i)

if(length(add_i) == 1)
    if(isempty(list_i))
        list_i = add_i;
    else
        list_i(end + 1) = add_i;
    end
else
    if(isempty(list_i))
        list_i = add_i;
    else
        for i = 1:length(add_i)
            list_i(end + 1) = add_i(i); %#ok<AGROW>
        end
    end
end

list_o = list_i;

end
