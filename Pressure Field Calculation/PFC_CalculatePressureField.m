function [Pressure_dxyz, Potential_dxyz, scan_axis, scan_list] = ...
    PFC_CalculatePressureField(...
    source_transducer_list, input_source,...                % tx setting
    reading_pos_x, reading_pos_y, reading_pos_z,...         % rx setting
    varargin)
% Calculate the Pressure Field using Rayleigh Sommerfeld Equation
% Subfunctions : this_VararginControl, this_VararginRearrange,
%                this_TriCenter, this_TriNV, this_TriArea, this_PFC_FrequencyDomain, 
%                this_PFC_FrequencyDomain_SinglePoint, this_ProgressBar

% [grid_varargin, mode_varargin, rx_varargin, one_varargin] = this_VararginControl(varargin{:});
[zero_varargin, one_varargin] = this_VararginControl(varargin{:});

one_varargin_order = one_varargin(1:2:end);
one_varargin_value = one_varargin(2:2:end);

%% Source Data
source_transducer = PFC_TransducerArray2Transducer(source_transducer_list);

Tri_Points = source_transducer.Tri_Points;
Tri_ConnectivityList = source_transducer.Tri_ConnectivityList;

source_coordinate = this_TriCenter(Tri_Points, Tri_ConnectivityList);
source_nv = this_TriNV(Tri_Points, Tri_ConnectivityList);
source_ds = this_TriArea(Tri_Points, Tri_ConnectivityList);

source_P = source_transducer.Pressure;

%% Input Data
input_type = input_source.Type;
reading_list = input_source.reading_list;

%% Reading Grid Data

[is_this_setting, Locb] = ismember({'meshgrid', 'ndgrid', 'nogrid'}, zero_varargin);

if(any(is_this_setting))
    grid_type = zero_varargin{max(Locb)};
else
    grid_type = 'meshgrid';
end

switch(grid_type)
    case 'meshgrid'
        [reading_pos_xxx, reading_pos_yyy, reading_pos_zzz] =...
            meshgrid(reading_pos_x, reading_pos_y, reading_pos_z);
        
    case 'ndgrid'
        [reading_pos_xxx, reading_pos_yyy, reading_pos_zzz] =...
            ndgrid(reading_pos_x, reading_pos_y, reading_pos_z);
        
    case 'nogrid'
        reading_pos_xxx = reading_pos_x;
        reading_pos_yyy = reading_pos_y;
        reading_pos_zzz = reading_pos_z;
end

%% Mode Data

[is_this_property, Locb] = ismember('split', one_varargin_order);

if(is_this_property)
    mode_type = 'split';
    split_num = one_varargin_value{Locb};
    
    source_coordinate_split_cnt = round(size(source_coordinate, 1)./split_num);
    scsc = source_coordinate_split_cnt;
else
    mode_type = 'full';
end

%% Rx Data

[is_this_property, Locb] = ismember('rx_transducer', one_varargin_order);

if(is_this_property)
    rx_info.is_rx_transducer = true;
    
    rx_transducer = one_varargin_value{Locb};
    
    rx_Tri_Points = rx_transducer.Tri_Points;
    rx_Tri_ConnectivityList = rx_transducer.Tri_ConnectivityList;
    rx_xyz = this_TriCenter(rx_Tri_Points, rx_Tri_ConnectivityList);
    rx_ds = this_TriArea(rx_Tri_Points, rx_Tri_ConnectivityList);
    assignin('base', 'rx_ds', rx_ds);
    rx_ds_ratio = rx_ds./sum(rx_ds, 'all');
    
    rx_info.rx_xyz = rx_xyz;
    rx_info.rx_ds_ratio = rx_ds_ratio;
else
    rx_info.is_rx_transducer = false;
    rx_info.rx_xyz = [0, 0, 0];
end

%% Progress Bar Data
[is_this_property, Locb] = ismember('update_percent', one_varargin_order);

if(is_this_property)
    pb_update_percent = one_varargin_value{Locb};
else
    pb_update_percent = 1;
end

%% Medium Data

[is_this_property, Locb] = ismember('medium', one_varargin_order);

if(is_this_property)
    Medium_name = one_varargin_value{Locb};
else
    Medium_name = 'air'; % Last version is set as water
end

medium_data = PFC_GetMedium('air'); % Last version is medium_name

%% GPU Check
try
    gpu_dev = gpuDevice;
catch
    gpu_dev = [];
end

switch(length(gpu_dev))
    case 0
    otherwise
        source_coordinate = gpuArray(source_coordinate);
        source_nv = gpuArray(source_nv);
        source_P = gpuArray(source_P);
        source_ds = gpuArray(source_ds);
        
        reading_pos_xxx = gpuArray(reading_pos_xxx);
        reading_pos_yyy = gpuArray(reading_pos_yyy);
        reading_pos_zzz = gpuArray(reading_pos_zzz);
        
        reading_list = gpuArray(reading_list);
end

%% Calculate Pressure

switch(input_type)
    case {'frequency', 'f'}
        last_axis = 'f';
        reading_f = reading_list;
        
        switch(mode_type)
            case 'full'
                this_ProgressBar('start', pb_update_percent,...
                    numel(reading_pos_xxx), size(rx_info.rx_xyz, 1), 1);
                
                [Pressure_dxyz, Potential_dxyz] = this_PFC_FrequencyDomain(...
                    source_coordinate, source_nv, source_ds, source_P,...
                    reading_pos_xxx, reading_pos_yyy, reading_pos_zzz, reading_f,...
                    medium_data, gpu_dev, rx_info, 1);
                
                this_ProgressBar('end');
                
            case 'split'
                Potential_dxyz = zeros([numel(reading_pos_xxx), length(reading_f)]);
                Pressure_dxyz = zeros([numel(reading_pos_xxx), length(reading_f)]);

                this_ProgressBar('start', pb_update_percent,...
                    numel(reading_pos_xxx), size(rx_info.rx_xyz, 1), split_num);
                
                for i = 1:split_num
                    if(i == split_num)
                        source_coordinate_i = source_coordinate((1 + scsc*(i - 1)):end, :);
                        source_nv_i = source_nv((1 + scsc*(i - 1)):end, :);
                        source_ds_i = source_ds((1 + scsc*(i - 1)):end, :);
                        source_P_i = source_P((1 + scsc*(i - 1)):end, :);
                    else
                        source_coordinate_i = source_coordinate((1:scsc) + scsc*(i - 1), :);
                        source_nv_i = source_nv((1:scsc) + scsc*(i - 1), :);
                        source_ds_i = source_ds((1:scsc) + scsc*(i - 1), :);
                        source_P_i = source_P((1:scsc) + scsc*(i - 1), :);
                    end
                    
                    [Pressure_dxyz_i, Potential_dxyz_i] = this_PFC_FrequencyDomain(...
                        source_coordinate_i, source_nv_i, source_ds_i, source_P_i,...
                        reading_pos_xxx, reading_pos_yyy, reading_pos_zzz, reading_f,...
                        medium_data, gpu_dev, rx_info, i);
                    
                    Potential_dxyz = Potential_dxyz + Potential_dxyz_i;
                    Pressure_dxyz = Pressure_dxyz + Pressure_dxyz_i;
                end
                
                this_ProgressBar('end');
        end
        
    case {'time', 't'}
        last_axis = 't';
        reading_time = reading_list;
        
        switch(mode_type)
            case 'full'
                this_ProgressBar('start', pb_update_percent,...
                    numel(reading_pos_xxx), size(rx_info.rx_xyz, 1), 1);
                
                [Pressure_dxyz, Potential_dxyz] = this_PFC_TimeDomain(...
                    source_coordinate, source_nv, source_ds, source_P, input_source,...
                    reading_pos_xxx, reading_pos_yyy, reading_pos_zzz, reading_time,...
                    medium_data, gpu_dev, rx_info, 1);
                
                this_ProgressBar('end');
                
            case 'split'
                Potential_dxyz = zeros([numel(reading_pos_xxx), length(reading_time)]);
                Pressure_dxyz = zeros([numel(reading_pos_xxx), length(reading_time)]);

                this_ProgressBar('start', pb_update_percent,...
                    numel(reading_pos_xxx), size(rx_info.rx_xyz, 1), split_num);
                
                for i = 1:split_num
                    if(i == split_num)
                        source_coordinate_i = source_coordinate((1 + scsc*(i - 1)):end, :);
                        source_nv_i = source_nv((1 + scsc*(i - 1)):end, :);
                        source_ds_i = source_ds((1 + scsc*(i - 1)):end, :);
                        source_P_i = source_P((1 + scsc*(i - 1)):end, :);
                    else
                        source_coordinate_i = source_coordinate((1:scsc) + scsc*(i - 1), :);
                        source_nv_i = source_nv((1:scsc) + scsc*(i - 1), :);
                        source_ds_i = source_ds((1:scsc) + scsc*(i - 1), :);
                        source_P_i = source_P((1:scsc) + scsc*(i - 1), :);
                    end
                    
                    [Pressure_dxyz_i, Potential_dxyz_i] = this_PFC_TimeDomain(...
                        source_coordinate_i, source_nv_i, source_ds_i, source_P_i, input_source,...
                        reading_pos_xxx, reading_pos_yyy, reading_pos_zzz, reading_time,...
                        medium_data, gpu_dev, rx_info, i);
                    
                    Potential_dxyz = Potential_dxyz + Potential_dxyz_i;
                    Pressure_dxyz = Pressure_dxyz + Pressure_dxyz_i;
                end
                
                this_ProgressBar('end');
        end
end




if(length(size(reading_pos_xxx)) == 3)
    Pressure_dxyz = reshape(Pressure_dxyz, [size(reading_pos_xxx), length(reading_list)]);
    Potential_dxyz = reshape(Potential_dxyz, [size(reading_pos_xxx), length(reading_list)]);
elseif(length(size(reading_pos_xxx)) == 2)
    Pressure_dxyz = reshape(Pressure_dxyz, [size(reading_pos_xxx), 1, length(reading_list)]);
    Potential_dxyz = reshape(Potential_dxyz, [size(reading_pos_xxx), 1, length(reading_list)]);
else
    
end








[is_this_setting] = ismember('scandata', zero_varargin);

if(is_this_setting)
    reading_list = gather(reading_list);
    
    switch(grid_type)
        case 'meshgrid'
            scan_axis_i = {'y', 'x', 'z', last_axis};
            scan_list_i = {reading_pos_y, reading_pos_x, reading_pos_z, reading_list};
            
        case 'ndgrid'
            scan_axis_i = {'x', 'y', 'z', last_axis};
            scan_list_i = {reading_pos_x, reading_pos_y, reading_pos_z, reading_list};
            
        case 'nogrid'
    end
    
    scan_axis = scan_axis_i;
    scan_list = scan_list_i;
    
%     for i = 1:length(scan_list_i)
%         scan_list_size(i) = length(scan_list_i{i}); %#ok<AGROW>
%     end
%     
%     scan_axis = [scan_axis_i(scan_list_size > 1), scan_axis_i(scan_list_size == 1)];
%     scan_list = [scan_list_i(scan_list_size > 1), scan_list_i(scan_list_size == 1)];
   

    
    scandata.scan_v = Pressure_dxyz;
    scandata.scan_axis = scan_axis;
    scandata.scan_list = scan_list;
    
    Pressure_dxyz = scandata;
end


end



%% Other Functions

function [zero_varargin, one_varargin] = this_VararginControl(varargin)

zero_varargin_list = {...
    'meshgrid', 'ndgrid', 'nogrid', 'scandata'};

one_varargin_list = {...
    'split', ...
    'rx_transducer',...
    'update_percent',...
    'medium'};

varargin_data = varargin;
varargin_list = {zero_varargin_list, one_varargin_list};
varargin_sub_number = [0, 1, 1];

varargin_data_rearranged = this_VararginRearrange(varargin_data, varargin_list, varargin_sub_number);

zero_varargin = varargin_data_rearranged{1};
one_varargin = varargin_data_rearranged{2};

end

function varargin_data_rearranged = this_VararginRearrange(varargin_data, varargin_list, varargin_sub_number)

varargin_data_rearranged = cell(size(varargin_list));
for i = 1:length(varargin_list)
    varargin_data_rearranged{i} = {};
end
varargin_data_rearranged{end + 1} = {};

v_count = 1;

while(v_count <= length(varargin_data))
    no_match = 1;
    for i = 1:length(varargin_list)
        try
            if(ismember(varargin_data(v_count), varargin_list{i}))
                varargin_data_rearranged{i}(end + 1) = varargin_data(v_count);
                for j = 1:varargin_sub_number(i)
                    v_count = v_count + 1;
                    varargin_data_rearranged{i}(end + 1) = varargin_data(v_count);
                end
                no_match = 0;
                break;
            end
        catch
            break;
        end
    end
    
    if(no_match)
        varargin_data_rearranged{end}(end + 1) = varargin_data(v_count);
    end
    
    v_count = v_count + 1;
end

end


function Center_Point = this_TriCenter(Points, ConnectivityList)
P = Points;
C = ConnectivityList;

P1 = P(C(:, 1), :);
P2 = P(C(:, 2), :);
P3 = P(C(:, 3), :);

Center_Point = (P1 + P2 + P3)./3;
end

function nv = this_TriNV(Points, ConnectivityList)
P = Points;
C = ConnectivityList;

P1 = P(C(:, 1), :);
P2 = P(C(:, 2), :);
P3 = P(C(:, 3), :);

v1 = P2 - P1;
v2 =  P3 - P2;

C = cross(v1, v2);
C_length = sqrt(sum(C.^2, 2));
C_length_matrix = C_length*ones(1, 3);
nv = C./C_length_matrix;
end

function A = this_TriArea(Points, ConnectivityList)
P = Points;
C = ConnectivityList;

P1 = P(C(:, 1), :);
P2 = P(C(:, 2), :);
P3 = P(C(:, 3), :);

a = sqrt(sum((P1 - P2).^2, 2));
b = sqrt(sum((P2 - P3).^2, 2));
c = sqrt(sum((P3 - P1).^2, 2));

s = (a + b + c)/2;

A = sqrt(s.*(s - a).*(s - b).*(s - c));
end



function [Pressure_dxyz, Potential_dxyz] = this_PFC_FrequencyDomain(...
    source_coordinate, source_nv, source_ds, source_P,...
    reading_pos_xxx, reading_pos_yyy, reading_pos_zzz, reading_f,...
    medium_data, gpu_dev, rx_info, split_num)

Potential_dxyz = zeros([numel(reading_pos_xxx), length(reading_f)]);
Pressure_dxyz = zeros([numel(reading_pos_xxx), length(reading_f)]);

is_rx_transducer = rx_info.is_rx_transducer;

if(is_rx_transducer)
    rx_xyz = rx_info.rx_xyz;
    rx_ds_ratio = rx_info.rx_ds_ratio;
    
    for ii = 1:numel(reading_pos_xxx)
        for jj = 1:size(rx_xyz, 1)
            reading_pos = [reading_pos_xxx(ii), reading_pos_yyy(ii), reading_pos_zzz(ii)]...
                + rx_xyz(jj, :);
            
            [Pressure_dxyz_i, Potential_dxyz_i] = this_PFC_FrequencyDomain_SinglePoint(...
                source_coordinate, source_nv, source_ds, source_P,...
                reading_pos, reading_f,...
                medium_data, gpu_dev);
            
            Potential_dxyz(ii, :) = Potential_dxyz(ii, :) + Potential_dxyz_i.*rx_ds_ratio(jj);
            Pressure_dxyz(ii, :) = Pressure_dxyz(ii, :) + Pressure_dxyz_i.*rx_ds_ratio(jj);
            
            this_ProgressBar('check', ii, jj, split_num);
        end
    end
else
    for ii = 1:numel(reading_pos_xxx)
        reading_pos = [reading_pos_xxx(ii), reading_pos_yyy(ii), reading_pos_zzz(ii)];
        
        [Pressure_dxyz_i, Potential_dxyz_i] = this_PFC_FrequencyDomain_SinglePoint(...
            source_coordinate, source_nv, source_ds, source_P,...
            reading_pos, reading_f,...
            medium_data, gpu_dev);
        
        Potential_dxyz(ii, :) = Potential_dxyz_i;
        Pressure_dxyz(ii, :) = Pressure_dxyz_i;
        
        this_ProgressBar('check', ii, 1, split_num);
    end
end
end

function [Pressure_dxyz, Potential_dxyz] = this_PFC_FrequencyDomain_SinglePoint(...
    source_coordinate, source_nv, source_ds, source_P,...
    reading_pos, reading_f,...
    medium_data, gpu_dev)

% size(source_coordinate)	[n, 3]
% size(source_nv)           [n, 3]
% size(source_ds)           [n, 1]
% size(source_P)            [n, 1]

% size(reading_pos)         [1, 3]
% size(reading_f)           [1, m]

% size(medium_info)         [1, 1]

%% Resize Input Data

switch(length(gpu_dev))
    case 0
        rho = medium_data.rho; % size is [1, 1]
        c = medium_data.c; % size is [1, 1]
        
        ones_n1 = ones([size(source_coordinate, 1), 1]); % size is [n, 1]
        ones_1m = ones([1, size(reading_f, 2)]); % size is [1, m]
    otherwise
        rho = gpuArray(medium_data.rho); % size is [1, 1]
        c = gpuArray(medium_data.c); % size is [1, 1]
        
        ones_n1 = gpuArray(ones([size(source_coordinate, 1), 1])); % size is [n, 1]
        ones_1m = gpuArray(ones([1, size(reading_f, 2)])); % size is [1, m]
end

reading_pos = ones_n1 * reading_pos; % size is [n, 3]
reading_f = ones_n1 * reading_f; % size is [n, m]

source_ds = source_ds * ones_1m; % size is [n, m]

r_vec = (reading_pos - source_coordinate); % size is [n, 3]
r = sqrt(sum(r_vec.^2, 2)); % size is [n, 1]
cos_theta = sum(r_vec.*source_nv, 2)./r; % size is [n, 1]

source_v = source_P./rho./c; % size is [n, 1]
source_v = source_v * ones_1m; % size is [n, m]

r = r * ones_1m; % size is [n, m]
cos_theta = cos_theta * ones_1m; % size is [n, m]

w = 2.*pi.*reading_f; % size is [n, m]
k = 2.*pi.*reading_f./(c); % size is [n, m]

phase_delay = exp(-1i.*k.*r); % size is [n, m]

velocity_potential_i = 1./(2.*pi).*cos_theta./r.*source_v.*phase_delay.*source_ds;
pressure_i = rho.*1i.*w.*velocity_potential_i;

velocity_potential = sum(velocity_potential_i, 1); % size is [1, m]
pressure = sum(pressure_i, 1); % size is [1, m]

Potential_dxyz = gather(velocity_potential);
Pressure_dxyz = gather(pressure);

end




function [Pressure_dxyz, Potential_dxyz] = this_PFC_TimeDomain(...
    source_coordinate, source_nv, source_ds, source_P, input_data,...
    reading_pos_xxx, reading_pos_yyy, reading_pos_zzz, reading_time,...
    medium_data, gpu_dev, rx_info, split_num)

rho = medium_data.rho; % size is [1, 1]
dt = reading_time(2) - reading_time(1);

Pressure_dxyz = zeros([numel(reading_pos_xxx), length(reading_time)]);
Potential_dxyz = zeros([numel(reading_pos_xxx), length(reading_time)]);

is_rx_transducer = rx_info.is_rx_transducer;

if(is_rx_transducer)
    rx_xyz = rx_info.rx_xyz;
    rx_ds_ratio = rx_info.rx_ds_ratio;
    
    for ii = 1:numel(reading_pos_xxx)
        for jj = 1:size(rx_xyz, 1)
            reading_pos = [reading_pos_xxx(ii), reading_pos_yyy(ii), reading_pos_zzz(ii)]...
                + rx_xyz(jj, :);
            
            [Potential_dxyz_i] = this_PFC_TimeDomain_SinglePoint(...
                source_coordinate, source_nv, source_ds, source_P, input_data,...
                reading_pos, reading_time,...
                medium_data, gpu_dev);
            
            Potential_dxyz(ii, :) = Potential_dxyz(ii, :) + Potential_dxyz_i.*rx_ds_ratio(jj);
            
            this_ProgressBar('check', ii, jj, split_num);
        end
    end
else
    for ii = 1:numel(reading_pos_xxx)
        reading_pos = [reading_pos_xxx(ii), reading_pos_yyy(ii), reading_pos_zzz(ii)];
        
        [Potential_dxyz_i] = this_PFC_TimeDomain_SinglePoint(...
            source_coordinate, source_nv, source_ds, source_P, input_data,...
            reading_pos, reading_time,...
            medium_data, gpu_dev);
        
        Potential_dxyz(ii, :) = Potential_dxyz_i;
        Pressure_dxyz(ii, :) = rho.*this_Gradient(Potential_dxyz_i, dt);
        
        this_ProgressBar('check', ii, 1, split_num);
    end
end
end

function [Potential_dxyz] = this_PFC_TimeDomain_SinglePoint(...
    source_coordinate, source_nv, source_ds, source_P, input_data,...
    reading_pos, reading_time,...
    medium_data, gpu_dev)

% size(source_coordinate)	[n, 3]
% size(source_nv)           [n, 3]
% size(source_ds)           [n, 1]
% size(source_P)            [n, 1]

% size(input_data)          [1, 1]

% size(reading_pos)         [1, 3]
% size(reading_time)        [1, m]

% size(medium_info)         [1, 1]

% size(gpu_dev)             [1, 1]

%% Resize Input Data
rho = medium_data.rho; % size is [1, 1]
c = medium_data.c; % size is [1, 1]

if(isfield(medium_data, 'alpha'))
    m_alpha = medium_data.alpha;
else
    m_alpha = 0;
end

switch(length(gpu_dev))
    case 0
        ones_n1 = ones([size(source_coordinate, 1), 1]); % size is [n, 1]
        ones_1m = ones([1, size(reading_time, 2)]); % size is [1, m]
    otherwise
        rho = gpuArray(rho); % size is [1, 1]
        c = gpuArray(c); % size is [1, 1]
        
        m_alpha = gpuArray(m_alpha);
        
        ones_n1 = gpuArray(ones([size(source_coordinate, 1), 1])); % size is [n, 1]
        ones_1m = gpuArray(ones([1, size(reading_time, 2)])); % size is [1, m]
end

reading_pos = ones_n1 * reading_pos; % size is [n, 3]
reading_time = ones_n1 * reading_time; % size is [n, m]

source_ds = source_ds * ones_1m; % size is [n, m]

r_vec = (reading_pos - source_coordinate); % size is [n, 3]
r = sqrt(sum(r_vec.^2, 2)); % size is [n, 1]
cos_theta = sum(r_vec.*source_nv, 2)./r; % size is [n, 1]

source_P_mag = abs(source_P) * ones_1m; % size is [n, m]

r = r * ones_1m; % size is [n, m]
cos_theta = cos_theta * ones_1m; % size is [n, m]
travel_time = r./c; % size is [n, m]

waveform_function = input_data.waveform_function;

travel_reading_time = reading_time - travel_time; % size is [n, m]

surface_pressure = waveform_function(travel_reading_time);
v_n_delay = surface_pressure./rho./c;

%% Medium Loss
med_loss = exp(-m_alpha.*r); % size is [n, m]

%%
vn_r = med_loss.*source_P_mag.*v_n_delay.*cos_theta./r.*source_ds; % size is [n, m]

Potential_dxyz = gather(1./(2.*pi).*sum(vn_r, 1)); % size is [1, m]

end




function this_ProgressBar(input_type, varargin)

global pb;
global pb_data;

switch(input_type)
    case 'start'
        try
            if(isempty(pb))
            else
                close(pb);
            end
        catch
        end
        
        pb = waitbar(0,'Start Setting');
        
        update_percent = varargin{1};
        
        input_data = [];
        
        for i = 2:length(varargin)
            input_data = [input_data, varargin{i}]; %#ok<AGROW>
        end

        data_size = input_data;
        if(length(data_size) == 1)
            data_size = [data_size, 1];
        end
        
        current_time = datetime('now');
        
        pb_data.update_percent = update_percent;
        pb_data.time_duration_matrix = zeros(data_size);
        pb_data.time_checked = false(data_size);
        
        pb_data.start_time = current_time;
        pb_data.total_num = prod(data_size);
        
        pb_data.current_time = current_time;
        pb_data.percent_floor = 0;
        pb_data.processed_time = 0;
        pb_data.processed_num = 0;
        
    case 'check'
        past_time = pb_data.current_time;
        past_percent_floor = pb_data.percent_floor;
        index_num = sub2ind(size(pb_data.time_duration_matrix), varargin{:});
        
        current_time = datetime('now');
        processed_time_i = seconds(current_time - past_time);
        
        pb_data.time_duration_matrix(index_num) = processed_time_i;
        pb_data.time_checked(index_num) = true;
        
        pb_data.current_time = current_time;
        pb_data.processed_time = pb_data.processed_time + processed_time_i;
        pb_data.processed_num = pb_data.processed_num + 1;
        
        average_time = pb_data.processed_time./pb_data.processed_num;
        process_left_num = pb_data.total_num - pb_data.processed_num;
        process_percentage = pb_data.processed_num/pb_data.total_num;
        time_to_go = average_time.*process_left_num;
        predicted_ending_time = current_time + duration(0,0,time_to_go);
        
        current_percent_floor = floor(100.*process_percentage/pb_data.update_percent);
        
        if(current_percent_floor > past_percent_floor)
            pb_data.percent_floor = current_percent_floor;
            
            try
                if(isempty(pb))
                    pb = waitbar(0,'Start Setting');
                end
                
                waitbar(process_percentage, pb, ...
                    {[num2str(100.*process_percentage, '%.1f'), '% (',...
                    num2str(pb_data.processed_num), '/', num2str(pb_data.total_num), ')',...
                    ' Expected Time Left : ', char(duration(0,0,time_to_go))], ...
                    ['Expected End Time : ', char(predicted_ending_time)]});
            catch
            end
        else
        end
    case 'end'
        start_time = pb_data.start_time;
        
        try
            if(isempty(pb))
            else
                close(pb);
            end
        catch
        end
        
        current_time = datetime('now');
        
        disp(['Started at ', char(start_time)]);
        disp(['Finished at ', char(current_time)]);
        disp(['Progressed time : ', char(current_time - start_time)]);
        
        clear global pb pb_data;
end

end


function F_g = this_Gradient(F, dt)
F_g = gradient(F, dt);
end