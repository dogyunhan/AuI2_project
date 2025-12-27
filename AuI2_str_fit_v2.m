clc; clearvars;

%% 1. 설정 및 데이터 로드
% -------------------------------------------------------------------------
% [User Input] 원자 번호 설정 (순서: [Ligand1, CentralAtom, Ligand2])
% 예: AuI2의 경우 I-Au-I 이므로 [53, 79, 53] (가운데가 중심 원자)
elem = [53, 79, 53]; 

% 상수 및 데이터 로드
run atom_consts.m % xfactor 변수 생성

% 데이터 경로 설정 (기존 유지)
base_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD";
dat_path = fullfile(base_path, "AuI2_30mM_0002", "PEPCed_dat.dat");
solv_path = fullfile(base_path, "heating_MeCN_0001", "merged_solv_dat.dat");
sads_path = fullfile(base_path, "AuI2_30mM_0002", "SADS_comps_3.dat");
std_sads_path = fullfile(base_path, "AuI2_30mM_0002", "std_SADS_comps_3.dat");

% 데이터 읽기
dat = readmatrix(dat_path);
% q = dat(:, 1);

% Solvent Heating Data Load
solv_data = readmatrix(solv_path);
q_w = solv_data(:, 1);
q_w = q_w(1<q_w & q_w<7);

heat_dat = solv_data(:, 2:4);
heat_dat = heat_dat(1<q_w & q_w<7, :); % Masking
[Uw, Sw, Vw] = svd(heat_dat, 'econ');
heat_dat = Uw(:, 1:3);

% SADS Data Load
SADS = readmatrix(sads_path);
std_SADS = readmatrix(std_sads_path);
q_sads = SADS(:, 1); % q vector for fitting

target_curve = SADS(:, 2); % 1st component (change index if needed)
target_std = std_SADS(:, 2);

%% 2. Scattering Factor 계산
% -------------------------------------------------------------------------
[f2, ff] = DHanfuncs.calc_scattering_factors(q_sads, elem, xfactor);

%% 3. Fitting 설정 (Config)
% -------------------------------------------------------------------------
cfg = struct();
cfg.q = q_sads;
cfg.target_dSq = target_curve;
cfg.target_Std = target_std;
cfg.f2 = f2;
cfg.ff = ff;
cfg.heat_dat = heat_dat;

% [General] Fitting Parameters: [r1, r2, theta]
% r1: Distance between Atom 1 and Atom 2 (Central)
% r2: Distance between Atom 3 and Atom 2 (Central)
cfg.para0 = [2.3, 2.3, 165];  % Initial guess
cfg.lb    = [2.0, 2.0, 150];    % Lower bounds
cfg.ub    = [3.0, 3.0, 180];  % Upper bounds
cfg.method = 'multistart';    

% Reference State (Ground State) Parameters for dSq calculation
% GS 구조: r1, r2, theta
% cfg.GS_params = [2.3542, 2.3542, 180]; % 학부생 ver

% [On the gold–ligand covalency in linear [AuX2]− complexes] ver
cfg.GS_params = [2.611,  2.611, 180];  % PBE(ADF)
% cfg.GS_params = [2.512, 2.512, 180];   % MP2
% cfg.GS_params = [2.540, 2.540, 180];   % SCS-MP2
% cfg.GS_params = [2.567, 2.567, 180];   % CCSD
% cfg.GS_params = [2.561, 2.561, 180];   % CCSD(T)

% [1,2,4-Triazole-derived carbene complexes of gold:
%  Characterization, solid-state aggregation and ligand disproportionation] ver
% cfg.GS_params = [2.63428, 2.63428, 180];   % 고체일때?


% Execute the fitting process and store the output
%% 4. Fitting 실행
% -------------------------------------------------------------------------
out = fit_Triatomic_Debye_PEPC(cfg);

%% 5. 결과 확인
% -------------------------------------------------------------------------
plot_data = struct();

plot_data(1).x = cfg.q; plot_data(1).y = cfg.target_dSq;
plot_data(1).color = 'red'; plot_data(1).label = 'Target SADS (PEPCed)';

plot_data(2).x = cfg.q; plot_data(2).y = out.fit_dSq_pepc;
plot_data(2).color = 'blue'; plot_data(2).label = 'Theory (PEPCed)';


title = ['Fit Result: r1=' num2str(out.paraFit(1), '%.5f') ...
       ', r2=' num2str(out.paraFit(2), '%.5f') ...
       ', \theta=' num2str(out.paraFit(3), '%.5f') char(176)];

DHanfuncs.custom_plot(plot_data, LineWidth=1.5, Title=title, XLim=[3, 7]);
disp('--- Fitting Parameters ---');
fprintf('r1 (Atom1-Central): %.5f\n', out.paraFit(1));
fprintf('r2 (Atom3-Central): %.5f\n', out.paraFit(2));
fprintf('theta (deg): %.5f\n', out.paraFit(3));
fprintf('chi sq: %.5f\n', out.chi2Fit);


%% ========================================================================
%  Functions
% =========================================================================
function sq = calc_Triatomic_Sq(q, r1, r2, theta, dwfs, f2, ff)
    % CALC_TRIATOMIC_SQ
    %   Generic geometry calculator for A-B-C type molecule
    %   B (Atom 2) is at Origin (0,0,0)
    %   A (Atom 1) is at distance r1
    %   C (Atom 3) is at distance r2, angle theta
    
        q = q(:);  % 열벡터로 만듦
        if (r1 <= 0) || (r2 <= 0), sq = ones(size(q))*realmax; return; end
        
        % Geometry Definition (Index 2 is Central Atom)
        % Atom 1 (Ligand 1)
        p1 = [-r1; 0; 0];
        % Atom 2 (Central)
        p2 = [0; 0; 0];
        % Atom 3 (Ligand 2)
        p3 = [-r2*cosd(theta); -r2*sind(theta); 0];
        
        % Calculate Distances
        % pdist order for 3 points: (1,2), (1,3), (2,3)
        % (1,2): r1 (Atom1-Central)
        % (1,3): Distance between Ligands
        % (2,3): r2 (Central-Atom3)
        r_vec = pdist([p1, p2, p3]'); 
        qr = q .* r_vec; % (N x 3)
        
        % Debye Equation
        % ff columns correspond to pdist pairs: (1,2), (1,3), (2,3)
        % Note: dwfs implementation simplified for now (can be expanded)
        
        sinc_qr = sin(qr) ./ qr;
        sinc_qr(qr==0) = 1;
        
        sq = f2 + 2 * sum(ff .* sinc_qr, 2);
    end

function out = fit_Triatomic_Debye_PEPC(cfg)
% FIT_TRIATOMIC_DEBYE_PEPC
%   일반화된 3원자 분자 피팅 함수
%   Fitting Variables: [r1, r2, theta]

    % Unpack
    q = cfg.q(:);
    target_dSq = cfg.target_dSq(:);
    target_Std = cfg.target_Std(:);
    f2 = cfg.f2;
    ff = cfg.ff;
    heat_dat = cfg.heat_dat;
    
    para0 = cfg.para0;
    lb = cfg.lb;
    ub = cfg.ub;
    
    % Reference (GS) Sq 계산
    GS_p = cfg.GS_params;
    GS_Sq = calc_Triatomic_Sq(q, GS_p(1), GS_p(2), GS_p(3), [0 0 0], f2, ff);
    
    % Objective Function
    objFun = @(p) chi2Obj(p, q, target_dSq, target_Std, f2, ff, GS_Sq, [0 0 0], heat_dat);
    
    % Optimization Options
    options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', true, 'MaxFunctionEvaluations', 1e5);
    problem = createOptimProblem('fmincon', 'x0', para0, 'objective', objFun, 'lb', lb, 'ub', ub, 'options', options);
    
    % Run Solver
    switch lower(cfg.method)
        case 'multistart'
            ms = MultiStart('UseParallel', true, 'Display', 'off');
            [paraFit, chi2Fit, ~, output] = run(ms, problem, 500); 
        case 'globalsearch'
            gs = GlobalSearch('NumTrialPoints', 500);
            [paraFit, chi2Fit, ~, output] = run(gs, problem);
        otherwise
            [paraFit, chi2Fit, ~, output] = fmincon(problem);
    end
    
    % Final Curves
    fit_Sq = calc_Triatomic_Sq(q, paraFit(1), paraFit(2), paraFit(3), [0 0 0], f2, ff);
    fit_dSq = fit_Sq - GS_Sq;
    [~, fit_dSq_pepc] = HKifuncs.pepc(fit_dSq, heat_dat); % Using external pepc function
    
    out.paraFit = paraFit;
    out.chi2Fit = chi2Fit;
    out.fit_dSq_pepc = fit_dSq_pepc;
    out.output = output;
end

function chi2 = chi2Obj(params, q, target_dSq, target_Std, f2, ff, GS_Sq, dwfs, heat_dat)
    r1 = params(1);
    r2 = params(2);
    theta = params(3);
    
    Sq = calc_Triatomic_Sq(q, r1, r2, theta, dwfs, f2, ff);
    % [~, Sq] = DHanfuncs.normalize_AUC(q, Sq, [4, 7]);
    % [~, GS_Sq] = DHanfuncs.normalize_AUC(q, GS_Sq, [4, 7]);
    
    theory_dSq = Sq - GS_Sq;
    [~, theory_dSq_pepc] = HKifuncs.pepc(theory_dSq, heat_dat);
    
    fit_idx = (q > 3);
    res = (target_dSq(fit_idx, :) - theory_dSq_pepc(fit_idx, :)) ./ target_Std(fit_idx, :);
    chi2 = sum(res.^2, 'omitnan');
end
