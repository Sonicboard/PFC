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
input_source.reading_list = 1e6; % [Hz]


%% Simulation Space
reading_pos_x = -5e-3 : 0.1e-3 : 5e-3;          % Simulation Area
reading_pos_y = -4e-3 : 0.1e-3 : 4e-3;          % Simulation Area
reading_pos_z = 5e-3 : 0.1e-3 : 60e-3;          % Simulation Area


%% Run Frequency Domain Simulation

% xz plane
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source,...                 % tx setting
reading_pos_x, 0, reading_pos_z);

Pressure_dxyz_i = abs(squeeze(Pressure_dxyz));

figure;imagesc(reading_pos_z.*1e3, reading_pos_x.*1e3, Pressure_dxyz_i);
set(gca,'DataAspectRatio',[1 1 1],'Layer','top');
colorbar;
xlabel('z [mm]');
ylabel('x [mm]');

% xy plane
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source,...                 % tx setting
reading_pos_x, reading_pos_y, 0.0166);

Pressure_dxyz_i = abs(squeeze(Pressure_dxyz));

figure;imagesc(reading_pos_x.*1e3, reading_pos_y.*1e3, Pressure_dxyz_i);
set(gca,'DataAspectRatio',[1 1 1],'Layer','top');
colorbar;
xlabel('x [mm]');
ylabel('y [mm]');

