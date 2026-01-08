clc; clearvars; close all;

%% ========================================================================
%  1. User Inputs & Configuration
% =========================================================================

% [System] 원자 번호 설정
elem_GS   = [53, 79, 53];  % Ground State: I-Au-I (Triatomic)
elem_prod = [79, 53];      % Product State: Au-I (Diatomic)
atom_I    = 53;            % Dissociated Atom: I

% [Path] 데이터 파일 경로
base_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD";
files = struct();
files.solv     = fullfile(base_path, "heating_MeCN_0001", "merged_solv_dat.dat");
files.sads     = fullfile(base_path, "AuI2_30mM_0002", "SADS_comps_4.dat"); 
files.sads_std = fullfile(base_path, "AuI2_30mM_0002", "std_SADS_comps_4.dat"); 

target_SADS = 1;
title = 'r_{Au-I} = %.4f A, r_{GS}1 = %.4f, r_{GS}2 = %.4f, theta = %.4f';

chi_red = false;

% [Fitting Parameters]
fit_range = [3.0, 7.0];    % q Fitting Range (A^-1)
init_pars = horzcat(2.3, [2.561 2.561 180]); 
lb        = horzcat(2.0, [2.45 2.45 180]);  % lower bound
ub        = horzcat(3.0, [2.65 2.65 180]);  % upper bound

% [External Script] 상수 로드
run atom_consts.m % xfactor 로드

%% ========================================================================
%  2. Data Loading & Preprocessing
% =========================================================================
fprintf('Loading and preprocessing data...\n');

% 2.1. Raw Data Load
raw_solv     = readmatrix(files.solv);
raw_dat     = readmatrix(files.sads);
raw_std = readmatrix(files.sads_std);

% 2.2. Define Master Mask (Slicing)
q_full = raw_dat(:, 1);
mask   = (q_full > 1) & (q_full < 7);

q_fit = q_full(mask);         % Fitting용 q 벡터
sads_comp = raw_dat(mask, target_SADS+1);    % idx=2이면 comp는 1
std_comp = raw_std(mask, target_SADS+1);     

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
[f2_prod, ff_prod] = DHanfuncs.calc_scattering_factors(q_fit, elem_prod, xfactor);

% Ground State (I-Au-I)
[f2_GS, ff_GS]     = DHanfuncs.calc_scattering_factors(q_fit, elem_GS, xfactor);

% Dissociated Atom (I)
[Sq_I, ~]          = DHanfuncs.calc_scattering_factors(q_fit, atom_I, xfactor);

%% ========================================================================
%  4. Fitting Configuration
% =========================================================================
cfg = struct();

cfg.fit_range = fit_range;
cfg.chi_red = chi_red;

% Data
cfg.q          = q_fit;
cfg.target_dSq = sads_comp;  % SADS comp는 이미 lsv 3개로 PEPC 되었음
cfg.target_Std = std_comp;
cfg.heat_dat   = heat_dat; % PEPC용 Basis

% Theory Factors
cfg.f2_prod = f2_prod;
cfg.ff_prod = ff_prod;
cfg.f2_GS   = f2_GS;
cfg.ff_GS   = ff_GS;
cfg.Sq_I    = Sq_I;

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
    r_prod = params(1);
    GS = [params(2), params(3), params(4)];  % r1, r2, theta
    
    % 1. Calculate Product State Sq (AuI + I)
    % Product is Diatomic (Au-I) + Monoatomic (I)
    Sq_prod = calc_Diatomic_Sq(cfg.q, r_prod, cfg.f2_prod, cfg.ff_prod);
    Sq_prod = Sq_prod + cfg.Sq_I; % Add dissociated Iodine atom
    
    % 2. Calculate Reference State Sq (AuI2)
    Sq_GS = calc_Triatomic_Sq(cfg.q, GS(1), GS(2), GS(3), cfg.f2_GS, cfg.ff_GS);
    
    % 3. Calculate Difference Spectrum (dSq)
    theory_dSq = Sq_prod - Sq_GS;
    
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