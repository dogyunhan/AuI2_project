clc; clearvars; close all;

%% ========================================================================
%  1. User Inputs & Configuration
% =========================================================================

% [System] 원자 번호 설정
elem_GS   = [53, 79, 53];  % Ground State: I-Au-I (Triatomic)
elem_iso = [79, 53, 53];      % Product State: Au-I-I (Isomer)
elem_AuI = [79, 53];      % Product State: Au-I
elem_I = 53;      % Product State: I

% [Path] 데이터 파일 경로
base_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD";
files = struct();
files.solv     = fullfile(base_path, "heating_MeCN_0001", "merged_solv_dat.dat");
files.PEPCed_dat = fullfile(base_path, "AuI2_30mM_0002", "PEPCed_dat.dat"); 
files.PEPCed_std = fullfile(base_path, "AuI2_30mM_0002", "PEPCed_std.dat"); 

times =  [0, 100e-12, 178e-12, 316e-12, 562e-12, ...
          1e-9, 1.78e-9, 3.16e-9, 5.62e-9, ...
          10e-9, 17.8e-9, 31.6e-9, 56.2e-9, ...
          100e-9, 178e-9, 316e-9, 562e-9, ...
          1e-6];
target_time = 100e-12;

title = ['r_{AuI} = %.4f / ' ...
    'r_{GS} = %.4f(2개), theta = %.4f / '];
chi_red = true;

% [Fitting Parameters]
fit_range = [3.0, 7.0];    % q Fitting Range (A^-1)
lb        = horzcat(2.4, [2.4 170]);  % lower bound
ub        = horzcat(2.9, [2.8 180]);  % upper bound
init_pars = lb;

% [External Script] 상수 로드
run atom_consts.m % xfactor 로드

%% ========================================================================
%  2. Data Loading & Preprocessing
% =========================================================================
fprintf('Loading and preprocessing data...\n');

% 2.1. Raw Data Load
raw_solv    = readmatrix(files.solv);
raw_dat     = readmatrix(files.PEPCed_dat);
raw_std     = readmatrix(files.PEPCed_std);

% 2.2. Define Master Mask (Slicing)
q_full = raw_dat(:, 1);
mask   = (q_full > 1) & (q_full < 7);

q_fit = q_full(mask);         % Fitting용 q 벡터
dat4time = raw_dat(mask, 3);    % idx=2이면 comp는 1
std4time = raw_std(mask, 3);     

% 2.3. Solvent Heating Data Processing
% 용매 데이터도 동일한 q grid를 갖는다고 가정하고 같은 mask 적용
q_solv = raw_solv(:, 1);
mask_solv = (q_solv > 1) & (q_solv < 7);
heat_dat = raw_solv(mask_solv, 2:end); 
[Uw, Sw, Vw] = svd(heat_dat, 'econ');
heat_dat = Uw(:, 1:3);
heat_dat = [heat_dat, ones(size(q_solv(mask_solv))), 1./q_solv(mask_solv)];

%% ========================================================================
%  3. Theory Calculation (Scattering Factors)
% =========================================================================
fprintf('Calculating scattering factors...\n');

% Product (Au-I)
[f2_AuI, ff_AuI] = DHanfuncs.calc_scattering_factors(q_fit, elem_AuI, xfactor);
[Sq_I, ~] = DHanfuncs.calc_scattering_factors(q_fit, elem_I, xfactor);

% Ground State (I-Au-I)
[f2_GS, ff_GS]     = DHanfuncs.calc_scattering_factors(q_fit, elem_GS, xfactor);

%% ========================================================================
%  4. Fitting Configuration
% =========================================================================
cfg = struct();

cfg.fit_range = fit_range;
cfg.chi_red = chi_red;

% Data
cfg.q          = q_fit;
cfg.target_dSq = dat4time;  % SADS comp는 이미 lsv 3개로 PEPC 되었음
cfg.target_Std = std4time;
cfg.heat_dat   = heat_dat; % PEPC용 Basis

% Theory Factors
cfg.Sq_I = Sq_I;
cfg.f2_AuI = f2_AuI;
cfg.ff_AuI = ff_AuI;
cfg.f2_GS   = f2_GS;
cfg.ff_GS   = ff_GS;

% Optimization Settings
cfg.x0     = init_pars;
cfg.lb     = lb;
cfg.ub     = ub;
cfg.method = 'multistart'; % 'multistart', 'globalsearch', or 'fmincon'

%% ========================================================================
%  5. Run Fitting
% =========================================================================
fprintf('Running optimization (%s)...\n', cfg.method);
out = run_structural_fitting(cfg);

%% ========================================================================
%  6. Visualization & Result
% =========================================================================
% Plotting
plot_data = struct();
plot_data(1).x = cfg.q; 
plot_data(1).y = cfg.target_dSq; 
plot_data(1).color = 'red'; 
plot_data(1).label = 'Experiment (PEPCed)';

plot_data(2).x = cfg.q; 
plot_data(2).y = out.fit_dSq;
plot_data(2).color = 'blue'; 
plot_data(2).label = 'Theory Fit';

plot_title = sprintf(title, out.params);

DHanfuncs.custom_plot(plot_data, LineWidth=1.5, Title=plot_title, XLim=[1 7]);

% Display Statistics
disp('========================================');
disp('           FITTING RESULTS              ');
disp('========================================');
fprintf('Fitting q-range: %.2f ~ %.2f\n', fit_range);
if chi_red
    fprintf('reduced Chi-squared value:   %.5f\n', out.chi2);
else
    fprintf('Chi-squared value:   %.5f\n', out.chi2);
end
disp('========================================');

%% ========================================================================
%  Local Functions
% =========================================================================

function out = run_structural_fitting(cfg)
    % Optimizer Setup
    objFun = @(p) objective_function(p, cfg);
    
    options = optimoptions('fmincon', 'Display', 'iter', ...
        'UseParallel', true, 'MaxFunctionEvaluations', 1e5);
    
    problem = createOptimProblem('fmincon', 'x0', cfg.x0, ...
        'objective', objFun, 'lb', cfg.lb, 'ub', cfg.ub, 'options', options);
    
    % Run Solver
    switch lower(cfg.method)
        case 'multistart'
            ms = MultiStart('UseParallel', true, 'Display', 'off');
            [p_opt, chi2, ~, output] = run(ms, problem, 500); 
        case 'globalsearch'
            gs = GlobalSearch('NumTrialPoints', 500);
            [p_opt, chi2, ~, output] = run(gs, problem);
        otherwise
            [p_opt, chi2, ~, output] = fmincon(problem);
    end
    
    % Generate Final Curve
    [~, fit_curve] = objective_function(p_opt, cfg);
    
    % Pack Output
    out.params  = p_opt; 
    out.chi2    = chi2;
    out.fit_dSq = fit_curve;
    out.output  = output;
end

function [chi2, theory_dSq_scaled] = objective_function(params, cfg)
    % Unpack
    r_AuI = params(1);
    GS = [params(2), params(2), params(3)];  % r1, r2, theta
    
    % 1. Calculate Product State Sq
    Sq_AuI = calc_Diatomic_Sq(cfg.q, r_AuI, cfg.f2_AuI, cfg.ff_AuI);

    % 2. Calculate Reference State Sq (AuI2)
    Sq_GS = calc_Triatomic_Sq(cfg.q, GS(1), GS(2), GS(3), cfg.f2_GS, cfg.ff_GS);
    
    % 3. Calculate Difference Spectrum (dSq)
    theory_dSq = (Sq_AuI + cfg.Sq_I) - Sq_GS;
    
    % 4. Apply PEPC & Scaling to match Experiment
    % (Orthogonalize against solvent heating)
    [~, theory_pepc] = HKifuncs.pepc(theory_dSq, cfg.heat_dat);
    [~, theory_dSq_scaled] = DHanfuncs.scaler(theory_pepc, cfg.target_dSq);

    % 5. Calculate Chi-square
    chi_mask = (cfg.q > cfg.fit_range(1)) & (cfg.q < cfg.fit_range(2));
    res = (cfg.target_dSq(chi_mask, :) - theory_dSq_scaled(chi_mask, :)) ./ cfg.target_Std(chi_mask, :);
    chi2 = sum(res.^2, 'omitnan');
    if cfg.chi_red
        chi2 = chi2 / (numel(cfg.q) - numel(params) - 1);
    end
end

function sq = calc_Diatomic_Sq(q, r, f2, ff)
    % Diatomic Debye Equation
    qr = q .* r;
    sinc_qr = sin(qr) ./ qr;
    sinc_qr(qr==0) = 1; % Handle division by zero
    sq = f2 + 2 * sum(ff .* sinc_qr, 2);
end

function sq = calc_Triatomic_Sq(q, r1, r2, theta, f2, ff)
    % Triatomic Debye Equation (General Geometry)
    % Atom 2 at Origin (0,0,0)
    
    if (r1 <= 0) || (r2 <= 0)
        error("Wrong Parameters! (should be larger than 0...")
    end
    
    % Coordinates
    p1 = [-r1; 0; 0];                              % Atom 1
    p2 = [0; 0; 0];                                % Atom 2 (Central)
    p3 = [-r2*cosd(theta); -r2*sind(theta); 0];    % Atom 3
    
    % Pairwise Distances: (1,2), (1,3), (2,3)
    r_vec = pdist([p1, p2, p3]'); 
    
    % Debye Sum
    qr = q .* r_vec;
    sinc_qr = sin(qr) ./ qr;
    sinc_qr(qr==0) = 1;
    
    sq = f2 + 2 * sum(ff .* sinc_qr, 2);
end