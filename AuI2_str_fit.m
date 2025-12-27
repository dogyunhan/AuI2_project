clc; clearvars;
run atom_consts.m
elem = [53, 79, 53]; % atom element number
% elem = [79, 53, 53]; % atom element number

dat = readmatrix(fullfile("\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuI2_30mM_0002", "PEPCed_dat.dat"));
q = dat(:, 1);

s = q/4/pi;  % sin(theta) / lambda (equal to q / 4Pi)

%% cross term 및 sqaure term 계산
tri_ff = zeros(3, numel(q));
for ii = 1:3
    for jj = 1:5
        tri_ff(ii, :) = tri_ff(ii, :) + (xfactor(elem(ii), jj).*exp(-xfactor(elem(ii), jj+6).*s.*s)');
    end
    tri_ff(ii, :) = tri_ff(ii, :) + xfactor(elem(ii), 6);
end
tri_ff = tri_ff';  % (3, q) -> (q, 3)
f2_raw = zeros(length(q), 3);
ff = zeros(numel(q), nchoosek(3, 2));  % 3C2개의 cross term 존재
cnt = 1;
for i = 1:3
    f2_raw(:, i) = tri_ff(:, i).^2;  % f2_raw = tri_ff * tri_ff
    for j = i+1:3
        ff(:, cnt) = tri_ff(:, i) .* tri_ff(:, j);  % cross term 계산
        cnt = cnt + 1;
    end
end
f2 = sum(f2_raw, 2);  % f2는 단순합 / ff는 sin(qr)/qr 곱해서 더해야 함

%%

solv_PATH = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\heating_MeCN_0001";
solv_data = readmatrix(fullfile(solv_PATH, "merged_solv_dat.dat"));

q_w = solv_data(:, 1);
heat_dat = solv_data(:, 2:4);

idx = 1<q_w & q_w<7;
heat_dat = heat_dat(idx, :);

%%
SADS_comp = 1;

SATS_PATH = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuI2_30mM_0002";
SADS = readmatrix(fullfile(SATS_PATH, "SADS_comps_3.dat"));  % pepced SADS (1~3)
std_SADS = readmatrix(fullfile(SATS_PATH, "std_SADS_comps_3.dat"));  % pepced SADS (1~3)

q = SADS(:, 1);
SADS = SADS(:, 2:end);
std_SADS = std_SADS(:, 2:end);

%% config 설정
cfg.q = q;
cfg.target_dSq = SADS(:, SADS_comp);
cfg.target_Std = std_SADS(:, SADS_comp);  % Assuming uniform weights
cfg.f2 = f2;
cfg.ff = ff;
cfg.heat_dat = heat_dat;
cfg.para0 = [1; 1; 180];  % Initial guess
cfg.lb = [0; 0; 0];  % Lower bounds
cfg.ub = [3; 3; 360];  % Upper bounds
cfg.method = 'multistart';  % fmincon, globalsearch, multistart 중 하나
% 
% cfg.GS_params = [2.3542 2.3542 180];
%%

out = fit_AuI2_Debye_PEPC(cfg);

%%
figure(1); clf; hold on;
plot(q, SADS(:, 1), 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'SADS');
plot(q, out.fit_dSq_pepc, 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'fit (PEPCed)');

disp(out.paraFit);

%%
function out = fit_AuI2_Debye_PEPC(cfg)
% FIT_AUI2_DEBYE_PEPC
%   AuI2- triatomic Debye model fitting:
%   GS_Sq fixed = d_AuI2_all(q, 2.611, 2.611, 180, 0,0,0, f2, ff)
%   Fit variables = [IAu, AuI, theta]
%   Objective = chi-square between PEPC(theory_dSq) and target_dSq
%
% Required fields in cfg:
%   cfg.q          (N×1) q grid
%   cfg.target_dSq (N×1) target curve (e.g. SADS(:,1))
%   cfg.target_Std (N×1) std (sigma) for chi-square weighting
%   cfg.f2         (N×1)
%   cfg.ff         (N×3) pair cross terms in pdist order (1,2)(1,3)(2,3)
%   cfg.heat_dat          2nd arg passed to HKifuncs.pepc (same as your script)
%
%   cfg.para0      (3×1) initial guess [IAu; AuI; theta_deg]
%   cfg.lb         (3×1) lower bound
%   cfg.ub         (3×1) upper bound
%   cfg.method     'fmincon' | 'globalsearch' | 'multistart'
%
% Optional:
%   cfg.dwfsFit    [dwfIAu dwfAuI dwfII] (default [0 0 0])
%   cfg.GS_params  [IAu AuI theta] (default [2.611 2.611 180])
%   cfg.options    optimoptions('fmincon',...)  (if not given, default created)
%   cfg.gsNumTrialPoints, cfg.gsNumStageOnePoints
%   cfg.nStartPoints (MultiStart)
%
% Output struct out:
%   out.paraFit, out.chi2Fit, out.GS_Sq, out.fit_Sq, out.fit_dSq, out.fit_dSq_pepc
%   out.solutions (solver outputs)

    % ---------- unpack & checks ----------
    q = cfg.q(:);
    target_dSq = cfg.target_dSq(:);
    target_Std = cfg.target_Std(:);

    f2 = cfg.f2(:);
    ff = cfg.ff;

    if size(ff,1) ~= numel(q) || size(ff,2) ~= 3
        error('cfg.ff must be N×3 where N = length(q).');
    end
    if numel(target_dSq) ~= numel(q) || numel(target_Std) ~= numel(q)
        error('cfg.target_dSq and cfg.target_Std must match length(cfg.q).');
    end

    para0 = cfg.para0(:);
    lb    = cfg.lb(:);
    ub    = cfg.ub(:);

    if numel(para0) ~= 3 || numel(lb) ~= 3 || numel(ub) ~= 3
        error('para0/lb/ub must be 3×1: [IAu; AuI; theta].');
    end

    method = lower(cfg.method);

    heat_dat = cfg.heat_dat;

    % ---------- DWF settings ----------
    dwfsFit = [0 0 0];
    if isfield(cfg,'dwfsFit') && ~isempty(cfg.dwfsFit)
        dwfsFit = cfg.dwfsFit;
    end

    % ---------- fixed GS ----------
    GS_params = [2.611 2.611 180];
    if isfield(cfg,'GS_params') && ~isempty(cfg.GS_params)
        GS_params = cfg.GS_params;
    end

    GS_Sq = d_AuI2_all(q, GS_params(1), GS_params(2), GS_params(3), 0, 0, 0, f2, ff);

    % ---------- objective ----------
    objFun = @(p) chi2Obj(p, q, target_dSq, target_Std, f2, ff, GS_Sq, dwfsFit, heat_dat);

    % ---------- fmincon options ----------
    if isfield(cfg,'options') && ~isempty(cfg.options)
        options = cfg.options;
    else
        options = optimoptions('fmincon', ...
            'Display','iter-detailed', ...      % 반복마다 상세 로그 출력
            'MaxFunctionEvaluations', 5e5, ...  % 목적함수 eval수
            'UseParallel', true, ...            % 병렬 계산
            'PlotFcn', @optimplotfval);         % 최적화 과정 중 목적함수 값 plot
    end

    problem = createOptimProblem('fmincon', ...
        'x0', para0, ...            % init parameter 설정
        'objective', objFun, ...    % 최소화할 목적함수(chi square)
        'lb', lb, ...               
        'ub', ub, ...
        'options', options);

    solutions = struct();

    % ---------- solve ----------
    switch method
        case 'fmincon'
            [paraFit, chi2Fit, exitflag, output] = fmincon(problem);  % exitflag: 종료 상태 알려줌 / output: 상세로그
            solutions.method   = 'fmincon';
            solutions.exitflag = exitflag;
            solutions.output   = output;

        case 'globalsearch'
            gsNumTrialPoints = 200;  
            gsNumStageOnePoints = 200;
            if isfield(cfg,'gsNumTrialPoints');      gsNumTrialPoints = cfg.gsNumTrialPoints; end
            if isfield(cfg,'gsNumStageOnePoints'); gsNumStageOnePoints = cfg.gsNumStageOnePoints; end

            gs = GlobalSearch('NumTrialPoints', gsNumTrialPoints, ...
                              'NumStageOnePoints', gsNumStageOnePoints);  % global search 객체 생성

            [paraFit, chi2Fit, exitflag, output, solset] = run(gs, problem);  % 여러 시작점에서 fmincon을 돌려봄(solset) => 최적값 반환(paraFit)
            solutions.method   = 'globalsearch';
            solutions.exitflag = exitflag;
            solutions.output   = output;
            solutions.solset   = solset;

        case 'multistart'
            nStartPoints = 50;  % 시작점 50개(default)
            if isfield(cfg,'nStartPoints'); nStartPoints = cfg.nStartPoints; end

            ms = MultiStart('UseParallel', true, 'Display', 'off');  % 객체 생성
            [paraFit, chi2Fit, exitflag, output, solset] = run(ms, problem, nStartPoints);
            solutions.method   = 'multistart';
            solutions.exitflag = exitflag;
            solutions.output   = output;
            solutions.solset   = solset;

        otherwise
            error("Unknown cfg.method: %s (use 'fmincon'|'globalsearch'|'multistart')", cfg.method);
    end

    % ---------- build final curves ----------
    % fitting한 para를 이용해 Sq를 계산
    fit_Sq = d_AuI2_all(q, paraFit(1), paraFit(2), paraFit(3), ...
                        dwfsFit(1), dwfsFit(2), dwfsFit(3), f2, ff);

    fit_dSq = fit_Sq - GS_Sq;

    [~, fit_dSq_pepc] = HKifuncs.pepc(fit_dSq, heat_dat);

    % ---------- outputs ----------
    out = struct();
    out.paraFit = paraFit;
    out.chi2Fit = chi2Fit;
    out.GS_Sq   = GS_Sq;
    out.fit_Sq  = fit_Sq;
    out.fit_dSq = fit_dSq;
    out.fit_dSq_pepc = fit_dSq_pepc;
    out.solutions = solutions;
end

%%
% ================= local functions =================

function chi2 = chi2Obj(params, q, target_dSq, target_Std, f2, ff, GS_Sq, dwfs, heat_dat)
    IAu   = params(1);
    AuI   = params(2);
    theta = params(3);

    Sq = d_AuI2_all(q, IAu, AuI, theta, ...
                    dwfs(1), dwfs(2), dwfs(3), f2, ff);

    theory_dSq = Sq - GS_Sq;
    [~, theory_dSq_pepc] = HKifuncs.pepc(theory_dSq, heat_dat);

    res  = (target_dSq - theory_dSq_pepc) ./ target_Std;
    chi2 = sum(res.^2);

    if ~isfinite(chi2)
        chi2 = realmax;
    end
end

%%
function sq = d_AuI2_all(q, IAu, AuI, theta, dwfIAu, dwfAuI, dwfII, f2, ff)
    q = q(:);
    sq = zeros(size(q));

    if (IAu <= 0) || (AuI <= 0)
        sq(:) = realmax;
        return
    end

    % atom1=I1, atom2=Au, atom3=I2
    pI1 = [-IAu; 0; 0];
    pAu = [0; 0; 0];
    pI2 = [-AuI*cosd(theta); -AuI*sind(theta); 0];

    r  = pdist([pI1, pAu, pI2]');      % (1,2)(1,3)(2,3)
    qr = q .* r;                       % N×3

    % pair order = (1,2)=IAu, (1,3)=II, (2,3)=AuI
    ff_dwf = exp(-q.^2 * [dwfIAu, dwfII, dwfAuI] / 2) .* ff;

    sinc_qr = sin(qr) ./ qr;
    sinc_qr(qr==0) = 1;                % 안전 처리

    sq = f2 + 2 * sum(ff_dwf .* sinc_qr, 2);
end
