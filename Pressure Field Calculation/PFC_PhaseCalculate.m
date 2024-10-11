function phase = PFC_PhaseCalculate(target, array, Lf)
%%
    % Phase_calculate - Computes clock counters and phase adjustments for an array of transducers.
    %
    % Syntax: [clk_counter, phase] = Phase_calculate(target, array, Lf, sos, fc, d)
    %
    % Inputs:
    %    target - Target position [x, y] in meters
    %    array - Size of the transducer array [rows, cols]
    %    Lf - Focal length in meters
    %    sos - Speed of sound in medium (m/s)
    %    fc - Center frequency of the transducer (Hz)
    %    d - Diameter between transducers in meters
    %
    % Outputs:
    %    clk_counter - Clock counters for each transducer
    %    phase - Phase adjustments for each transducer

    % Basic Setting
    sos = 343;
    fc = 40000;
    d = 0.01;
    r = d / 2;

    % Calculate wavelength
    lambda = sos / fc;

    % Preallocate arrays for time and clock counters
    time = zeros(array(1), array(2));
    clk_counter = zeros(array(1), array(2));

    % Calculate time of flight and clock counter for each transducer
    for i = 1:array(1)
        for j = 1:array(2)
            % Calculate the distance from each transducer to the target
            R1 = norm([(i-1) * d + r, (j-1) * d + r] - target);
            time(i, j) = sqrt(R1^2 + Lf^2) / sos;
            phase(i, j) = 2*pi*fc*time(i,j);
        end
    end

