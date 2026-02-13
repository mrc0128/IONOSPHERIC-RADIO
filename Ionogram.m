setenv('DIR_MODELS_REF_DAT', '/Users/micahchannell/pharlap/dat');

clear; close all; clc;

%% Paramter setup
% Locations
auburn_lat = 32.6;   % Reciever
auburn_lon = -85.5;

chesap_lat = 36.8;  % transmitter 1
chesap_lon = -76.3;

corpus_lat = 27.8;   % transmitter 2
corpus_lon = -97.4;

% Micah Day
year = 2003;
month = 6;
day = 24;
hour = 24;  % UTC
Time = [year, month, day, hour, 0];

% Frequency range
freq_start = 1;      % MHz
freq_end = 10;       % MHz
freq_step = 0.1;     % MHz

%fprintf('Receiver: Auburn, AL (%.2f, %.2f)\n', auburn_lat, auburn_lon);
%fprintf('Transmitter 1: Chesapeake, VA (%.2f, %.2f)\n', chesap_lat, chesap_lon);
%fprintf('Transmitter 2: Corpus, Texas (%.2f, %.2f)\n', corpus_lat, corpus_lon);
%fprintf('Date and Time: %d-%02d-%02d %02d:00 \n', year, month, day, hour);

% calculate range from auburn to chesap

orig_lat = chesap_lat;
orig_lon = chesap_lon;
[range_ches, bearing_ches] = latlon2raz(auburn_lat, auburn_lon, orig_lat, ...
    orig_lon, 'wgs84');

fprintf('\nChesapeake to Auburn:\n');
fprintf('    Range: %.2f km\n', range_ches/1000);
fprintf('    Bearing: %.2f degrees\n', bearing_ches);

% calculate range from auburn to corpus 
orig_lat = corpus_lat;
orig_lon = corpus_lon;
[range_corpus, bearing_corpus] = latlon2raz(auburn_lat, auburn_lon, orig_lat, ...
    orig_lon, 'wgs84');

fprintf('\nCorpus to Auburn:\n');
fprintf('    Range: %.2f km\n', range_corpus/1000);
fprintf('    Bearing: %.2f degrees\n', bearing_corpus);



R12 = 63.6; % sunspot 
max_range = 5000; % num range step (increased 10%)
num_range = 201; 
range_inc = max_range/(num_range-1);
start_height = 0;
height_inc = 2;
num_heights = 500; % num of height steps
kp = 0; % kp index 
doppler_flag = 1; % doppler parameter 


% generate ionogram for chesapeake transmitter
tic
[iono_pf_grid_ches, iono_pf_grid_5_ches, collision_freq_ches, irreg_ches, iono_te_grid_ches] = gen_iono_grid_2d(chesap_lat, chesap_lon, R12, Time, ...
    bearing_ches, max_range, num_range, range_inc, ...
    start_height, height_inc, num_heights, kp, doppler_flag);
toc

fprintf('Grid Size: %d heights x %d ranges\n', num_heights, num_range);

iono_en_grid_ches = iono_pf_grid_ches.^2 / 80.6164e-6;
iono_en_grid_5_ches = iono_pf_grid_ches.^2 / 80.6164e-6;

% ray trace at 5 MHz

fprintf('\nSingle ray trace for 5MHz\n');

% ray param

test_freq = 5.0; %MHz
test_elev = 10; % degrees (first guess)
tol = [1e-7 0.01 25]; % ODE, target, iterations 
nhops = 1; %hops
ir_flag = 0; % irregularities flag 

% trace ray 

tic

[ray_data, ray_path_data, ray_state_vec] = raytrace_2d(chesap_lat, chesap_lon, test_elev, bearing_ches, test_freq, nhops, tol, ir_flag, iono_en_grid_ches, iono_en_grid_5_ches, collision_freq_ches, start_height, height_inc, range_inc, irreg_ches);
toc

fprintf('final ground range: %.2f km\n', ray_data.ground_range);

fprintf('Find elevation angle');
target_range = range_ches/1000; %km
freq = 5; %MHz
elevs = 0.5:0.5:60; %elevation grid

% find good rays at all elevation
freqs = freq * ones(size(elevs));
[ray_data_fan, ~, ~] = raytrace_2d(chesap_lat, chesap_lon, elevs, bearing_ches, freqs, nhops, tol, ir_flag, iono_en_grid_ches, iono_en_grid_5_ches, collision_freq_ches, start_height, height_inc, range_inc, irreg_ches);


gnd_ranges = zeros(size(elevs));
for idx = 1:length(elevs)
    if ray_data_fan(idx).ray_label == 1
        gnd_ranges(idx) = ray_data_fan(idx).ground_range;
    else 
        gnd_ranges(idx) = NaN;
    end
end

valid_idx = find(~isnan(gnd_ranges));
fprintf(' %d good rays out of %d total\n', length(valid_idx), length(elevs));
%
fprintf('min range: %.2f km at %.1f degrees\n', min(gnd_ranges(valid_idx)), min(elevs(valid_idx(1))));
fprintf('max range: %.2f km at %.1f degrees\n', max(gnd_ranges(valid_idx)), max(elevs(valid_idx(end)))); 


% find which ray his target range
[min_diff, closest_idx] = min(abs(gnd_ranges(valid_idx)-target_range));

fprintf('closest ray is %.2f km away from target\n', min_diff);
fprintf('best elevation: %.1f degrees\n', elevs(valid_idx(closest_idx)));


% best ray idx

best_idx = valid_idx(closest_idx);

% calc height for ray

group_range = ray_data_fan(best_idx).group_range;
ground_range = ray_data_fan(best_idx).ground_range;
virtual_ht = ray_data_fan(best_idx).apogee;

fprintf('n\5 MHz: \n');
fprintf('   Ground Range: %.2f km\n', ground_range);
fprintf('   Group Range: %.2f km\n', group_range);
fprintf('   Virtual Height: %.2f km\n', virtual_ht);
fprintf('   Best Ray Index: %d\n', best_idx);
fprintf('   Elevation Angle: %.1f degrees\n', elevs(best_idx));

% frequency array
freqs_array = freq_start:freq_step:freq_end;
num_freqs = length(freqs_array);

% more arrays
virtual_ht_ches = NaN(num_freqs, 1);
elevations_ches = NaN(num_freqs, 1);

fprintf('Full Ionogram For Chesapeake');
fprintf('testing freq from %.1f to %.1f MHz\n', freq_start, freq_end);

% loop each frequency
for freq_idx = 1:num_freqs
    current_freq = freqs_array(freq_idx);

% trace rays at all elevations for freq

freq_vec = current_freq * ones(size(elevs));
[ray_data_fan, ~, ~] = raytrace_2d(chesap_lat, chesap_lon, elevs, ...
    bearing_ches, freq_vec, nhops, tol, ir_flag);

% extract ground ranges for rays

gnd_ranges = NaN(size(elevs));
for idx = 1:length(elevs)
    if ray_data_fan(idx).ray_label == 1
        gnd_ranges(idx) = ray_data_fan(idx).ground_range;
    end
end

% find valid rays
valid_idx = find(~isnan(gnd_ranges));

if isempty(valid_idx)
    fprintf('you got no rays lol');
else
    [min_diff, closest_idx] = min(abs(gnd_ranges(valid_idx)-target_range));
fprintf('The closest ray is %.2f km away from the target\n', min_diff);
fprintf('Best elevation: %.1f degrees\n', elevs(valid_idx(closest_idx)));

% store the results
    best_idx = valid_idx(closest_idx);
    elevations_ches(freq_idx) = elevs(best_idx);
    virtual_ht_ches(freq_idx) = ray_data_fan(best_idx).apogee;
    ground_ranges_ches(freq_idx) = ray_data_fan(best_idx).ground_range;
    group_ranges_ches(freq_idx) = ray_data_fan(best_idx).group_range;
end

end


% store array for ground ranges
%ground_ranges_ches = NaN(num_freqs, 1);
%group_ranges_ches = NaN(num_freqs, 1);

figure('Position', [100, 100, 1000, 800]);


% plot virtual height vs frequency 
subplot(2,2,1);
plot(freqs_array, virtual_ht_ches, 'b.-', 'LineWidth',2, 'MarkerSize',8);
grid on;
xlabel('Frequency (MHz)', 'FontSize', 11);
ylabel('Virtual Height (km)', 'FontSize', 11);
title(sprintf('Ionogram: Chesapeark to Auburn\nRange: %.0f km, Bearing: %.1f', target_range, bearing_ches), 'FontSize', 12);
xlim([freq_start freq_end]);

% plot elevation angle vs frequency

subplot(2,2,2);
plot(freqs_array, elevations_ches, 'r.-', 'LineWidth',2, 'MarkerSize',8);
grid on;
xlabel('Frequency (MHz', 'FontSize',11);
ylabel('Elevation Angle (°)', 'FontSize', 11);
title('Required Elevation Angle', 'FontSize', 12);
xlim([freq_start freq_end]);
ylim([0 60]);


% plot ground range vs frequency
subplot(2, 2, 3);
plot(freqs_array, ground_ranges_ches, 'g.-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
yline(target_range, 'k--', 'LineWidth', 1.5, 'DisplayName','Target');
grid on;
xlabel('Frequency (MHz)', 'FontSize', 11);
ylabel('Ground Range (km)', 'FontSize', 11);
title(' Acheived Ground Range', 'FontSize', 12);
xlim([freq_start freq_end]);
legend('Location','best');

% ground range vs frequency
subplot(2, 2, 4);
plot(freqs_array, group_ranges_ches, 'm.-', 'LineWidth', 2, 'MarkerSize',8 );
grid on;
xlabel('Frequency (MHz)', 'FontSize', 11);
ylabel('Group Range (km)', 'FontSize', 11);
title('Group Range (Virtual Path Length)', 'FontSize', 12);
xlim([freq_start freq_end]);