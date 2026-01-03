clc; clearvars; close all;

%% ========================================================================
%  1. User Inputs & Configuration
% =========================================================================

% [System] 원자 번호 설정
elem_bent = [53, 79, 53];
atom_Au = 79; 
atom_I = 53; 

% [Path] 데이터 파일 경로
base_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD";
files = struct();
files.solv     = fullfile(base_path, "heating_MeCN_0001", "merged_solv_dat.dat");
files.dads     = fullfile(base_path, "AuI2_30mM_0002", "DADS_comps_4.dat"); 
files.dads_std = fullfile(base_path, "AuI2_30mM_0002", "std_DADS_comps_4.dat"); 

target_DADS = 2;
title = 'r_{bent} = %.4f, %.4f, theta = %.4f';

% [Fitting Parameters]
fit_range = [1.0, 7.0];    % q Fitting Range (A^-1)
init_pars = horzcat([2.5 2.5 150]); 
lb        = horzcat([2.0 2.0 90]);  % lower bound
ub        = horzcat([3.0 3.0 180]);  % upper bound

% [External Script] 상수 로드
run atom_consts.m % xfactor 로드

%% ========================================================================
%  2. Data Loading & Preprocessing
% =========================================================================
fprintf('Loading and preprocessing data...\n');

% 2.1. Raw Data Load
raw_solv    = readmatrix(files.solv);
raw_dat     = readmatrix(files.dads);
raw_std     = readmatrix(files.dads_std);

% 2.2. Define Master Mask (Slicing)
q_full = raw_dat(:, 1);
mask   = (q_full > 1) & (q_full < 7);

q_fit = q_full(mask);         % Fitting용 q 벡터
dads_comp = (-1)*raw_dat(mask, target_DADS+1);    % idx=2이면 comp는 1
std_comp = raw_std(mask, target_DADS+1);     

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

[f2_bent, ff_bent] = DHanfuncs.calc_scattering_factors(q_fit, elem_bent, xfactor);

[Sq_Au, ~]  = DHanfuncs.calc_scattering_factors(q_fit, atom_Au, xfactor);
[Sq_I, ~]  = DHanfuncs.calc_scattering_factors(q_fit, atom_I, xfactor);


%% ========================================================================
%  4. Fitting Configuration
% =========================================================================
cfg = struct();

cfg.fit_range = fit_range;

% Data
cfg.q          = q_fit;
cfg.target_dSq = dads_comp;  % SADS comp는 이미 lsv 3개로 PEPC 되었음
cfg.target_Std = std_comp;
cfg.heat_dat   = heat_dat; % PEPC용 Basis

cfg.f2_bent = f2_bent;
cfg.ff_bent = ff_bent;

cfg.Sq_Au = Sq_Au;
cfg.Sq_I = Sq_I;

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
fprintf('Chi-squared value:   %.5f\n', out.chi2);
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
    BENT = [params(1), params(2), params(3)];  % r1, r2, theta
    
    Sq_bent = calc_Triatomic_Sq(cfg.q, BENT(1), BENT(2), BENT(3), cfg.f2_bent, cfg.ff_bent);
    
    % 3. Calculate Difference Spectrum (dSq)
    theory_dSq =(cfg.Sq_Au + 2*cfg.Sq_I) - Sq_bent;
    
    % 4. Apply PEPC & Scaling to match Experiment
    % (Orthogonalize against solvent heating)
    [~, theory_pepc] = HKifuncs.pepc(theory_dSq, cfg.heat_dat);
    [~, theory_dSq_scaled] = DHanfuncs.scaler(theory_pepc, cfg.target_dSq);

    % 5. Calculate Chi-square
    chi_mask = (cfg.q > cfg.fit_range(1)) & (cfg.q < cfg.fit_range(2));
    res = (cfg.target_dSq(chi_mask, :) - theory_dSq_scaled(chi_mask, :)) ./ cfg.target_Std(chi_mask, :);
    chi2 = sum(res.^2, 'omitnan');
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