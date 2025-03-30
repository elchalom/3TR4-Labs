load('lab4_num_expt3.mat'); % Load data (xt, yt, tt)

% Compute cross-correlation
[Rxy, lags] = xcorr(yt, xt); % Cross-correlation of y(t) and x(t)

% Find the delay T
[~, idx] = max(abs(Rxy)); % Peak of cross-correlation
tau = lags(idx); % Lag in samples
dt = tt(2) - tt(1); % Time step (sampling interval)
T_estimate = tau * dt; % Convert lag to seconds

% --- Plot 1: Input Signal x(t) ---
figure;
plot(tt, xt, 'b', 'LineWidth', 1.5);
title('Input Signal x(t)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- Plot 2: Output Signal y(t) ---
figure;
plot(tt, yt, 'r', 'LineWidth', 1.5);
title('Output Signal y(t) (Noisy)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- Plot 3: Cross-Correlation Rxy(τ) ---
figure;
plot(lags * dt, Rxy, 'k', 'LineWidth', 1.5); % Convert lags to time units
hold on;
plot(T_estimate, Rxy(idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Mark peak
title('Cross-Correlation R_{xy}(\tau)');
xlabel('Lag \tau (s)');
ylabel('Correlation');
legend('Cross-Correlation', 'Peak (Estimated T)');
grid on;

fprintf('Estimated delay T: %.4f seconds\n', T_estimate);