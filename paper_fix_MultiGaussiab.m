%% ========================================================================
% ULTIMATE MULTI-GAUSSIAN DECONVOLUTION & SCHERRER CRYSTALLITE SIZE ALGORITHM
% (UPDATED: FIXED FILL POLYGON & ALLOWED NEGATIVE VALUES)
% ========================================================================
clc; clear; close all;

fprintf('=================================================================\n');
fprintf('  MULTI-GAUSSIAN XRD DECONVOLUTION & SCHERRER ANALYSIS TOOL      \n');
fprintf('=================================================================\n\n');

%% 1. DATA IMPORT
% -------------------------------------------------------------------------
filename = '/MATLAB Drive/New Folder/Kak Nunu data/0.25.xlsx'; % Sesuaikan path Excel Anda

if exist(filename, 'file')
    data = readmatrix(filename);
    twoTheta = data(:, 1);
    I_raw    = data(:, 2);
    fprintf('[1/8] Data XRD berhasil diimpor dari file: %s\n', filename);
else
    fprintf('[1/8] File tidak ditemukan. Menggunakan Data DUMMY (Sagu+ZnO)...\n');
    twoTheta = (10:0.05:70)';
    I_amorf_dummy = 800 * exp(-((twoTheta - 20)/6).^2) + 400 * exp(-((twoTheta - 30)/8).^2);
    I_kristal_dummy = 1200*exp(-((twoTheta-31.8)/0.3).^2) + 1500*exp(-((twoTheta-34.4)/0.35).^2) + 900*exp(-((twoTheta-36.3)/0.3).^2);
    bg_dummy = 150 + 2*twoTheta;
    noise_dummy = 15 * randn(size(twoTheta));
    I_raw = I_amorf_dummy + I_kristal_dummy + bg_dummy + noise_dummy;
end

N = length(twoTheta);

%% 2. BACKGROUND CORRECTION (ALS) & SMOOTHING
% -------------------------------------------------------------------------
fprintf('[2/8] Koreksi Baseline Latar Belakang (ALS) & Smoothing...\n');
lambda = 1e5; p_als = 0.002;
I_bg = als_baseline(I_raw, lambda, p_als);

% Mengurangkan baseline dari data mentah
I_nobg = I_raw - I_bg; 
% PERBAIKAN: Baris I_nobg(I_nobg < 0) = 0; DIHAPUS agar nilai minus tetap terhitung

% Savitzky-Golay Filter untuk menghaluskan data yang sudah dikurangi background
I_smooth = sgolayfilt(I_nobg, 3, 21); 

%% 3. PEAK DETECTION & INITIAL PARAMETER ESTIMATION
% -------------------------------------------------------------------------
fprintf('[3/8] Deteksi Puncak & Estimasi Parameter Awal...\n');
MinPeakProm = 0.03 * max(I_smooth);
[pks, locs, widths, ~] = findpeaks(I_smooth, twoTheta, ...
    'MinPeakProminence', MinPeakProm, 'MinPeakDistance', 0.5);

num_peaks = length(pks);
p0 = zeros(1, 3 * num_peaks);
lb = zeros(1, 3 * num_peaks);
ub = zeros(1, 3 * num_peaks);

for i = 1:num_peaks
    idx = (i-1)*3;
    A_init = pks(i);
    mu_init = locs(i);
    sigma_init = widths(i) / 2.35482; % FWHM ke Sigma
    
    p0(idx+1) = A_init;       p0(idx+2) = mu_init;         p0(idx+3) = sigma_init;
    lb(idx+1) = -Inf;         lb(idx+2) = mu_init - 1.0;   lb(idx+3) = 0.05; % lb tinggi puncak diubah agar mentolerir lekukan
    ub(idx+1) = A_init * 1.5; ub(idx+2) = mu_init + 1.0;   ub(idx+3) = 15.0;
end

%% 4. MULTI-GAUSSIAN NON-LINEAR FITTING
% -------------------------------------------------------------------------
fprintf('[4/8] Optimasi Multi-Gaussian (Levenberg-Marquardt)...\n');
options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off', 'MaxFunctionEvaluations', 20000, 'MaxIterations', 1000);

p_opt = lsqcurvefit(@multi_gaussian_model, p0, twoTheta, I_smooth, lb, ub, options);

%% 5. AREA CALCULATION, CLASSIFICATION, & SCHERRER EQUATION
% -------------------------------------------------------------------------
fprintf('[5/8] Kalkulasi FWHM Analitik, Klasifikasi, dan Persamaan Scherrer...\n');
FWHM_threshold = 2.5; % Batas amorf vs kristalin (derajat)
K_scherrer = 0.9;     
lambda_xrd = 0.15406; % nm (Cu-Kalpha)

A_c_total = 0;
A_a_total = 0;
peak_results = struct();
peak_shapes = zeros(N, num_peaks);

for i = 1:num_peaks
    idx = (i-1)*3;
    A_i     = p_opt(idx+1);
    mu_i    = p_opt(idx+2);
    sigma_i = p_opt(idx+3);
    
    % Profil dan Luas Area Analitik
    peak_shapes(:, i) = A_i * exp(-((twoTheta - mu_i).^2) / (2 * sigma_i^2));
    fwhm_i = 2.35482 * sigma_i;
    area_i = A_i * sigma_i * sqrt(2 * pi);
    
    % Klasifikasi dan Ukuran Kristalit
    if fwhm_i <= FWHM_threshold
        classification = "Crystalline";
        A_c_total = A_c_total + area_i;
        
        theta_rad = deg2rad(mu_i / 2);
        beta_rad  = deg2rad(fwhm_i);
        D_size_nm = (K_scherrer * lambda_xrd) / (beta_rad * cos(theta_rad));
    else
        classification = "Amorphous";
        A_a_total = A_a_total + area_i;
        D_size_nm = NaN; 
    end
    
    % Simpan ke struktur
    peak_results(i).Peak_No = i;
    peak_results(i).Center_2Theta = mu_i;
    peak_results(i).Height = A_i;
    peak_results(i).FWHM_deg = fwhm_i;
    peak_results(i).Area = area_i;
    peak_results(i).Classification = classification;
    peak_results(i).Crystallite_Size_nm = D_size_nm;
end

I_fit_total = sum(peak_shapes, 2);
X_c = (A_c_total / (A_c_total + A_a_total)) * 100;
X_a = 100 - X_c;

%% 6. MENAMPILKAN HASIL DI COMMAND WINDOW
% -------------------------------------------------------------------------
fprintf('\n=================================================================\n');
fprintf('   HASIL KRISTALINITAS & UKURAN KRISTAL (HYBRID MULTI-GAUSSIAN)  \n');
fprintf('=================================================================\n');
fprintf('Persentase Kristalin (Xc) : %.5f %%\n', X_c);
fprintf('Persentase Amorf     (Xa) : %.5f %%\n', X_a);
fprintf('=================================================================\n');
fprintf('%-6s | %-10s | %-10s | %-12s | %-10s | %-11s\n', ...
    'Peak', '2-Theta', 'FWHM(deg)', 'Class', 'Area', 'Size(nm)');
fprintf('-----------------------------------------------------------------\n');
for i = 1:num_peaks
    size_str = sprintf('%.2f', peak_results(i).Crystallite_Size_nm);
    if isnan(peak_results(i).Crystallite_Size_nm), size_str = '-'; end
    fprintf('%-6d | %-10.3f | %-10.3f | %-12s | %-10.1f | %-11s\n', ...
        i, peak_results(i).Center_2Theta, peak_results(i).FWHM_deg, ...
        peak_results(i).Classification, peak_results(i).Area, size_str);
end
fprintf('=================================================================\n\n');

%% 7. VISUALISASI SIAP PUBLIKASI (2 SUBPLOT)
% -------------------------------------------------------------------------
fprintf('[6/8] Membuat Plot Visualisasi...\n');
figure('Color', 'w', 'Name', 'XRD Multi-Gaussian & Scherrer Analysis', 'Position', [100, 50, 950, 750]);

% Subplot 1: Dekonvolusi Multi-Gaussian
subplot(2,1,1); hold on;
plot(twoTheta, I_raw, 'k.', 'MarkerSize', 5, 'DisplayName', 'Raw Data');
plot(twoTheta, I_bg, 'g--', 'LineWidth', 1.2, 'DisplayName', 'ALS Background');
plot(twoTheta, I_fit_total + I_bg, 'r-', 'LineWidth', 2, 'DisplayName', 'Multi-Gaussian Fit');

for i = 1:num_peaks
    if peak_results(i).Classification == "Crystalline"
        col = [0 0.447 0.741]; % Biru
    else
        col = [0.85 0.325 0.098]; % Merah Bata
    end
    
    % PERBAIKAN LOGIKA FILL (Area diarsir dengan rapi mengikuti kontur)
    fill_x = [twoTheta; flipud(twoTheta)];
    fill_y = [(peak_shapes(:, i) + I_bg); flipud(I_bg)];
    
    fill(fill_x, fill_y, col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

title(sprintf('Dekonvolusi Multi-Gaussian (Xc = %.2f%%, Xa = %.2f%%)', X_c, X_a), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Intensitas (a.u.)', 'FontWeight', 'bold');
xlim([min(twoTheta) max(twoTheta)]); 

% PERBAIKAN Y-LIMIT: Mengakomodasi nilai minus secara dinamis
min_y_val = min([min(I_raw), min(I_bg)]);
max_y_val = max(I_raw);
ylim([min_y_val - abs(min_y_val*0.05), max_y_val + abs(max_y_val*0.1)]); 
grid on; grid minor;
legend('Location', 'northeast', 'Box', 'off');
set(gca, 'LineWidth', 1.2, 'TickDir', 'in');

% Subplot 2: Isolasi Puncak Kristalin & Label Scherrer
subplot(2,1,2); hold on;
I_cryst_only = zeros(N, 1);
for i = 1:num_peaks
    if peak_results(i).Classification == "Crystalline"
        I_cryst_only = I_cryst_only + peak_shapes(:, i);
    end
end
plot(twoTheta, I_cryst_only, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Fase Kristalin Diekstrak');

for i = 1:num_peaks
    if peak_results(i).Classification == "Crystalline"
        x0 = peak_results(i).Center_2Theta;
        h  = peak_results(i).Height;
        d  = peak_results(i).Crystallite_Size_nm;
        plot(x0, h, 'r*', 'MarkerSize', 6, 'HandleVisibility', 'off');
        text(x0, h + max(I_cryst_only)*0.08, sprintf('%.1f^\\circ\n(%.1fnm)', x0, d), ...
            'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    end
end
title('Estimasi Ukuran Kristalit (D) - Persamaan Scherrer', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Sudut 2\theta (degree)', 'FontWeight', 'bold'); ylabel('Intensitas (a.u.)', 'FontWeight', 'bold');
xlim([min(twoTheta) max(twoTheta)]); 

min_cryst_y = min(I_cryst_only);
ylim([min_cryst_y - abs(min_cryst_y*0.1), max(I_cryst_only)*1.2]);
grid on; grid minor;
set(gca, 'LineWidth', 1.2, 'TickDir', 'in');

%% 8. EKSPOR KE EXCEL
% -------------------------------------------------------------------------
fprintf('[7/8] Mengekspor Hasil ke Excel...\n');
output_file = 'Hasil_Lengkap_MultiGaussian_Scherrer.xlsx';
try
    T_peaks = struct2table(peak_results);
    writetable(T_peaks, output_file, 'Sheet', 'Peak_Details');
    
    T_sum = table(A_c_total, A_a_total, X_c, X_a, lambda_xrd, K_scherrer, ...
        'VariableNames', {'Area_Crystalline', 'Area_Amorphous', 'Crystallinity_Percent', ...
        'Amorphous_Percent', 'Xray_Wavelength_nm', 'Scherrer_K_Factor'});
    writetable(T_sum, output_file, 'Sheet', 'Summary');
    fprintf('[8/8] Analisis Selesai. Data diekspor ke: %s\n', output_file);
catch ME
    fprintf('Ekspor Excel Gagal: %s\n', ME.message);
end

%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================
function Y = multi_gaussian_model(p, x)
    Y = zeros(size(x));
    num_peaks = length(p) / 3;
    for i = 1:num_peaks
        idx = (i-1)*3;
        A = p(idx+1); mu = p(idx+2); sigma = p(idx+3);
        Y = Y + A * exp(-((x - mu).^2) / (2 * sigma^2));
    end
end

function z = als_baseline(y, lambda, p)
    m = length(y);
    D = diff(speye(m), 2);
    w = ones(m, 1);
    for iter = 1:10
        W = spdiags(w, 0, m, m);
        C = chol(W + lambda * (D' * D));
        z = C \ (C' \ (w .* y));
        w = p * (y > z) + (1 - p) * (y <= z);
    end
end