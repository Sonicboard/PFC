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
input_source_f = [];

input_source_f.Type = 'frequency';
input_source_f.reading_list = 1e6; % [Hz]


%% Time Domain
input_source_time = [];

input_source_time.Type = 'time';
input_source_time.waveform_function = @this_waveform_function;

T = 1./transducer_info.Frequency;
dt = T/100;
reading_time = (50*T):dt:(53*T);

input_source_time.reading_list = reading_time; % [s]


%% Simulation Space
reading_pos_x = -5e-3 : 0.1e-3 : 5e-3;          % Simulation Area
reading_pos_y = -4e-3 : 0.1e-3 : 4e-3;          % Simulation Area
reading_pos_z = 5e-3 : 0.1e-3 : 60e-3;          % Simulation Area


%% Run Frequency Domain Simulation

% z axis
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source_f,...                 % tx setting
0, 0, reading_pos_z);

Pressure_dxyz_frequency = abs(squeeze(Pressure_dxyz));

figure;plot(reading_pos_z.*1e3, Pressure_dxyz_frequency);

% z axis
[Pressure_dxyz] = ...
PFC_CalculatePressureField(...
transducer_object, input_source_time,...                 % tx setting
0, 0, reading_pos_z);

Pressure_dxyz_time = squeeze(Pressure_dxyz);

hold on;plot(reading_pos_z.*1e3, max(Pressure_dxyz_time, [], 2));



%%
function waveform_o = this_waveform_function(time_list)
% this_function should always have 'time_list' input and 'waveform_o'
% output

waveform_o = zeros(size(time_list));

f = 1e6;
T = 1/f;

time_start = 0;
time_end = 80.*T;

time_on = (time_list >= time_start) & (time_list <= time_end);

waveform_o(time_on) = sin(2.*pi.*f.*time_list(time_on));

end


