%% Transducer Design
transducer_info = [];

transducer_info.Frequency = 1e6;
transducer_info.Pressure = 1;
transducer_info.Phase = 0;
transducer_info.Type = 'circle';
transducer_info.Radius = 5e-3;

transducer_object = PFC_Make3DTransducer(transducer_info);

PFC_DrawTransducer3D(transducer_object, '-');


%% Frequency Domain
input_source = [];

input_source.Type = 'frequency';
reading_f = (0.1:0.01:5).*1e6; % [Hz]
input_source.reading_list = reading_f;


%% Simulation Space
reading_pos_x = -5e-3 : 0.1e-3 : 5e-3;          % Simulation Area
reading_pos_y = -4e-3 : 0.1e-3 : 4e-3;          % Simulation Area
reading_pos_z = 5e-3 : 0.1e-3 : 60e-3;          % Simulation Area


%% Run Frequency Domain Simulation

% xz plane
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source,...                 % tx setting
0, 0, reading_pos_z);

Pressure_dxyz_i = abs(squeeze(Pressure_dxyz));
scan1_v = Pressure_dxyz_i;

figure;imagesc(reading_f.*1e-6, reading_pos_z.*1e3, Pressure_dxyz_i);
xlabel('f [MHz]');
ylabel('z [mm]');


