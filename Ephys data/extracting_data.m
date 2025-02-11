%Reading of nachbac data

% wt_fly_1_CS = abfload('Fly_1_CS_0000.abf');
% wt_fly_2_CS = abfload('Fly_2_CS_0000.abf');
% wt_fly_3_CS = abfload('Fly_3_CS_0000.abf');
% wt_fly_4_CS = abfload('Fly_4_CS_0000.abf');
% wt_fly_5_CS = abfload('Fly_5_CS_0000.abf');
% wt_fly_6_CS = abfload('Fly_6_CS_0000.abf');
% wt_fly_7_CS = abfload('Fly_7_CS_0000.abf');
% wt_fly_8_CS = abfload('Fly_8_CS_0000.abf');
% wt_fly_9_CS = abfload('Fly_9_CS_0000.abf');
% wt_fly_10_CS = abfload('Fly_10_CS_0000.abf');
% wt_fly_11_CS = abfload('Fly_11_CS_0000.abf');
% wt_fly_12_CS = abfload('Fly_12_CS_0000.abf');

[num_points, num_channels, num_reps] = size(wt_fly_12_CS);

% Generate time points
timepoints = 0:2:(2*(num_points-1));

% Create a new figure
figure;

% Define the channel to plot (in this case, the first channel)
channel = 1;

% Loop through each repetition
for rep = 1:num_reps
    % Extract the trace for the current channel and repetition
    trace = squeeze(wt_fly_12_CS(:, channel, rep));
    
    % Plot the trace with timepoints on the x-axis
    plot(timepoints/100, trace);
    hold on; % Keep the current plot to overlay the next trace
end

% Add labels and title
xlabel('Time (s)');
ylabel('Value');
title('Fly 1 - Channel 1');

% Adjust layout
hold off;
xlabel('Time (ms)');
ylabel('Voltage (mV)');


combined_data = zeros(num_points, num_reps + 1);

% Assign timepoints to the first column of the matrix
combined_data(:, 1) = timepoints(:);

% Loop through each repetition and add the traces to the matrix
for rep = 1:num_reps
    % Extract the trace for the current channel and repetition
    trace = squeeze(wt_fly_12_CS(:, channel, rep));
    
    % Add the trace to the next column in the matrix
    combined_data(:, rep + 1) = trace;
end

% Define the filename for saving
filename = 'fly_12_wt_control.csv';

% Save the combined data to a CSV file
writematrix(combined_data, filename);

