%% Transducer Design
% tx transducer
transducer_info = [];

transducer_info.Frequency = 1e6;
transducer_info.Pressure = 1;
transducer_info.Phase = 0;
transducer_info.Type = 'circle';
transducer_info.Radius = 5e-3;

transducer_object = PFC_Make3DTransducer(transducer_info);

PFC_DrawTransducer3D(transducer_object, '-');

% rx transducer
transducer_info = [];

transducer_info.Frequency = 1e6;
transducer_info.Pressure = 1;
transducer_info.Phase = 0;
transducer_info.Type = 'circle';
transducer_info.Radius = 0.5e-3;
transducer_info.rotate = [180, 0, 0];

transducer_object_rx = PFC_Make3DTransducer(transducer_info);

PFC_DrawTransducer3D(transducer_object_rx, '-');


%% Frequency Domain
input_source = [];

input_source.Type = 'frequency';
input_source.reading_list = 1e6; % [Hz]


%% Simulation Space
reading_pos_x = -5e-3 : 0.1e-3 : 5e-3;          % Simulation Area
reading_pos_y = -4e-3 : 0.1e-3 : 4e-3;          % Simulation Area
reading_pos_z = 5e-3 : 0.1e-3 : 60e-3;          % Simulation Area


%% Run Frequency Domain Simulation

% z axis
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source,...                 % tx setting
0, 0, reading_pos_z);

Pressure_dxyz_frequency = abs(squeeze(Pressure_dxyz));

figure;plot(reading_pos_z.*1e3, Pressure_dxyz_frequency);

% z axis
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source,...                 % tx setting
0, 0, reading_pos_z, 'rx_transducer', transducer_object_rx);

Pressure_dxyz_frequency = abs(squeeze(Pressure_dxyz));

hold on;plot(reading_pos_z.*1e3, Pressure_dxyz_frequency);

