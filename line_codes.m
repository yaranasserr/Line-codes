% (1) Generate random bits of zeros and ones.
% (2) Modulate this same vector using the different types of line codes.
%NON RETURN TO ZERO/ NON RETURN TO ZERO INVERTED /RETURN TO ZERO/ AMI/ MANCHESTER/MULTILEVEL TARNSMISSION 3
% (3) Plot a sample of all previous line code modulation and plot them under each other on the
% same figure using subplot, ex subplot (6,1), (make sure that all the types have the same periodic time Ts).
% (4) Find the power spectrum density of each code, and plot them in the same figure using subplot as previous.
% Generate a random sequence of 0s and 1s
bits = randi([0 1], 1, 10);

% NON RETURN TO ZERO
modulated_bits = 2*bits - 1; % Convert 0s to -1s and 1s to +1s
figure(1);subplot(6, 1, 1);
stairs(modulated_bits,'LineWidth', 2);
ylim([-1.5 1.5]);
title('Non return to zero');
xlabel('Time');
ylabel('Amplitude');
% figure(3);

figure(2);plot_power_spectral_density(modulated_bits,1,'Power Spectral Density of NON RETURN TO ZERO')
%NON RETURN TO ZERO INVERTED 
%watch this https://www.youtube.com/watch?v=Kxndom8GaUQ
% Initialize the encoded signal with the first bit
lastBit = 1;
modulated_bits = zeros(1, length(bits));

% Perform NRZI encoding
for i = 1:length(bits)
    if bits(i) == 0
        modulated_bits(i) = lastBit;
    else
        modulated_bits(i) = -lastBit;
        lastBit = -lastBit;
    end
end


% Plot the encoded signal
figure(1);subplot(6, 1, 2);
stairs(modulated_bits,'LineWidth', 2);
ylim([-1.5 1.5]);
title('Non return to zero inverted');
xlabel('Time');
ylabel('Amplitude');
figure(2);plot_power_spectral_density(modulated_bits,2,'Power Spectral Density of Non return to zero inverted')
%Return to Zero
modulated_bits = zeros(1, 2*length(bits));

for i = 1:length(bits)
    if bits(i) == 0
        modulated_bits(2*i-1:2*i) = [-1, 0];
    else
        modulated_bits(2*i-1:2*i) = [1, 0];
    end
end

% Plot the modulated bits
figure(1);subplot(6, 1, 3);
t = 0:0.5:(length(modulated_bits)/2 - 0.5);
stairs(t, modulated_bits, 'LineWidth', 2);
ylim([-1.5 1.5]);
title('Return to Zero');
xlabel('Time');
ylabel('Amplitude');
figure(2);plot_power_spectral_density(modulated_bits,3,'Power Spectral Density of Return to Zero')


% AMI
% Initialize the AMI encoder
modulated_bits = zeros(1, length(bits));
last_voltage = 1;

% Perform AMI encoding
for i = 1:length(bits)
    if bits(i) == 0
        modulated_bits(i) = 0;
    else
        if last_voltage == 1
            modulated_bits(i) = 1;
            last_voltage = -1;
        else
            modulated_bits(i) = -1;
            last_voltage = 1;
        end
    end
end

% Plot the modulated bits
figure(1);subplot(6, 1, 4);
t = 0:1:length(modulated_bits)-1;
stairs(t, modulated_bits, 'LineWidth', 2);
ylim([-1.5 1.5]);
title('Modulated Bits using AMI');
xlabel('Time');
ylabel('Amplitude');
figure(2);plot_power_spectral_density(modulated_bits,4,'Power Spectral Density of AMI encoding')
%Manchester

% Initialize the Manchester encoder
modulated_bits = zeros(1, length(bits)*2);
for i = 1:length(bits)
    if bits(i) == 0
        modulated_bits((i-1)*2+1:i*2) = [1 -1];
    else
        modulated_bits((i-1)*2+1:i*2) = [-1 1];
    end
end

figure(1);subplot(6, 1, 5);
t = 0:0.5:length(modulated_bits)/2-0.5;
stairs(t, modulated_bits, 'LineWidth', 2);
ylim([-1.5 1.5]);
title('Modulated Bits using Manchester Coding');
xlabel('Time');
ylabel('Amplitude');
figure(2);plot_power_spectral_density(modulated_bits,5,'Power Spectral Density of Manchester Coding')

%MLT3
signal_levels = zeros(1, length(bits));
signal_levels(bits == 1) = 1;
signal_levels(bits == 0) = 0;
signal_levels = 2*signal_levels - 1;

% Apply multilevel transmission 3 modulation
% Apply multilevel transmission 3 modulation
modulated_bits = repelem(signal_levels, 3);

% Create a time vector for the modulated bit sequence
t = linspace(0, length(modulated_bits)/3, length(modulated_bits));

% Plot the modulated bit sequence
figure(1);subplot(6, 1, 6);
stairs(t, modulated_bits, 'LineWidth', 2);
ylim([-1.5 1.5]);
title('Multilevel Transmission 3 Modulation');
xlabel('Time');
ylabel('Signal Level');

figure(2);plot_power_spectral_density(modulated_bits,6,'Power Spectral Density of Multilevel Transmission 3 Modulation')
%Power spectral density 
function plot_power_spectral_density(modulated_bits,n,label)
% Compute the power spectral density of the modulated bits using periodogram
fs = length(modulated_bits); % sampling frequency
[Pxx, f] = periodogram(modulated_bits, [], [], fs);
%[Pxx, f] = pwelch(modulated_bits, [], [], [], fs);
% Plot the power spectral density
subplot(6, 1, n);
plot(f, Pxx, 'LineWidth', 2);
xlim([0 fs/2]);
title(label);
xlabel('Frequency');
ylabel('Power');
    end
