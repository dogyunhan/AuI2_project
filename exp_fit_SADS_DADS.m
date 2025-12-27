clc; clearvars;

default_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuI2_30mM_0002";
% default_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuBr2_150mM_0001";

data_merge_all = readmatrix(fullfile(default_path, "PEPCed_dat.dat"));
std_merge_all = readmatrix(fullfile(default_path, "PEPCed_std.dat"));

q = data_merge_all(:, 1);
data_merge_all = data_merge_all(:, 2:end);
std_merge_all = std_merge_all(:, 2:end);

% ps 단위로 이루어짐
tds_merge = [0, ...
          100, 178, 316, 562, ...
          1000, 1780, 3160, 5620, ...
          10000, 17800, 31600, 56200, ...
          100000, 178000, 316000, 562000, ...
          1000000];
%% =========================================
%  Std analysis of Vm (time components)
%  - comps_p: 사용할 SVD 성분 인덱스
% =========================================
[Up, Sp, Vp] = svd(data_merge_all, 'econ');
comps_p = [1 2 3];
[std_p, conc_p] = HKifuncs.fit_V_std(comps_p, data_merge_all, std_merge_all, Up, Sp, Vp, 10, tds_merge); %#ok<NASGU>

%% (옵션) 초기 t < 1 ps의 가중치 조정
% % 1 ps에서는 std를 크게 주어서 exp fit을 느슨하게 함
% tt_lim = 1;
% rel_weight = 100;
% temp_tind = find(tds_merge < tt_lim);
% std_p(temp_tind, :) = std_p(temp_tind, :) * rel_weight;


%% =========================================
%  Global exponential fitting (IRF estimation)
%  - 공유 IRF와 공유 시간상수(τ)로 여러 Vm 성분 동시 피팅
% =========================================
num_exps  = 5;
% num_exps  = 4;
num_sines = 0;

low_b  = [ 0 100];
up_b   = [ 0 100];
init_par = [0 100];

% ps 단위로

lims_tl = [0.02 800 10000 1e5 1e10];     % τ 하한
lims_tu = [0.02 1500 30000 7e5 1e10];  % τ 상한

% lims_tl = [0.02 1000 60000 1e10];     % τ 하한
% lims_tu = [0.02 1500 70000 1e10];  % τ 상한

lims_sl = [0.1 0.1 0.1 0.1];
lims_su = [2   2   2   2  ];

lims_rl = [0.02 0.02 0.02 0.02];
lims_ru = [5    5    5    5   ];

lims_dl = [0.02 0.02 0.02 0.02];
lims_du = [100  100  100  100 ];

num_Vs = numel(comps_p);

% 각 지수항의 성분별 계수(±500) 및 τ 초기값/경계 추가
for e = 1:num_exps-1
    for v = 1:num_Vs
        low_b   = [low_b, -500]; %#ok<AGROW>
        up_b    = [up_b,   500]; %#ok<AGROW>
        init_par = [init_par, 0]; %#ok<AGROW>
    end
    low_b   = [low_b, lims_tl(e)]; %#ok<AGROW>
    up_b    = [up_b,  lims_tu(e)]; %#ok<AGROW>
    init_par = [init_par, 0.5*(lims_tl(e) + lims_tu(e))]; %#ok<AGROW>
end
low_b   = [low_b, lims_tl(num_exps)];
up_b    = [up_b,  lims_tu(num_exps)];
init_par = [init_par, 0.5*(lims_tl(num_exps) + lims_tu(num_exps))];

% (진동항 없음) num_sines = 0 이므로 아래 루프는 실행 안됨
for s = 1:num_sines %#ok<UNRCH>
    for v = 1:num_Vs
        low_b   = [low_b, -500]; %#ok<AGROW>
        up_b    = [up_b,   500]; %#ok<AGROW>
        init_par = [init_par, 0]; %#ok<AGROW>
    end
    low_b   = [low_b, lims_sl(s), lims_rl(s), lims_dl(s), 0]; %#ok<AGROW>
    up_b    = [up_b,  lims_su(s), lims_ru(s), lims_du(s), 2.0*pi]; %#ok<AGROW>
    init_par = [init_par, ...
        0.5*(lims_sl(s)+lims_su(s)), ...
        0.5*(lims_rl(s)+lims_ru(s)), ...
        0.5*(lims_dl(s)+lims_du(s)), pi]; %#ok<AGROW>
end

% 예시 초기값 (최근 실험값으로 갱신 권장)
% init_par = [ ...
%   0.274641403531613, 0.295033740285038, ...
%  -0.0548809277857306, 0.0240288037464530, -0.101803156762615, ...
%   0.0200000000000000, -0.0302421380917050, 0.00703133884976953, ...
%   0.205394514179481, 3.34835979310419, ...
%  -0.0832155885294696, -0.658951589481646, -0.264386359287487, ...
%  1246.86789554659, 100000];

[theory_pl, theory_pl_exp, xl, errsl, fvall, hessianl] = HKifuncs.exp_fit5( ...
    Vp(:, comps_p), std_p(:, comps_p), tds_merge, ...
    num_exps, num_sines, ...
    low_b, up_b, ...
    lims_tl, lims_tu, lims_sl, lims_su, lims_dl, lims_du, ...
    init_par, [1 1e6], 20, false); %#ok<ASGLU>

% 결과 백업
theory_pl_bk     = theory_pl.o';
theory_pl_exp_bk = theory_pl_exp.o';

% 핵심 파라미터 인덱스 구성(t0, IRF, 각 τ들)
ind = [1 2];
for e = 1:num_exps-1
    ind = [ind, ind(end) + num_Vs + 1]; %#ok<AGROW>
end
ind = [ind, ind(end) + 1];
for s = 1:num_sines %#ok<UNRCH>
    ind = [ind, (ind(end)+num_Vs+1) : (ind(end)+num_Vs+4)]; %#ok<AGROW>
end

clc;
fprintf('%.6e ', xl(ind));
fprintf("\n");
fprintf('%.6e ', errsl(ind));
% disp(errsl(ind))
% disp(fvall)


%% SADS 뽑아내기
time_constants = xl([10 14 18 19]);
% time_constants = xl([8 11 14 15]);
% time_constants = xl([8 11 12]);


% [SAC, std_SAC, theory_profile, profile_SAC, std_profile_SAC] = DHanfuncs.KCA_for_AuI2(data_merge_all, std_merge_all, q, tds_merge, -1, time_constants, tds_merge, 111, [1, 1e6]);
[SAC, std_SAC, theory_profile, profile_SAC, std_profile_SAC] = HKifuncs.KCA(data_merge_all, std_merge_all, q, tds_merge, -1, time_constants, tds_merge, 111, [1, 1e6]);
% [DADS, std_DADS, theory_profile_d, profile_SAC_d, std_profile_SAC_d] = HKifuncs.KCA_DADS(data_merge_all, std_merge_all, q, tds_merge, -1, time_constants, tds_merge, 200, [1 1e6]);
%% SADS comp 고른 뒤 저장
SADS_comps = 3;

save = true;
if save
    writematrix([q SAC(:, 1:SADS_comps)], fullfile(default_path, sprintf("SADS_comps_%d.dat", SADS_comps)));
    writematrix([q std_SAC(:, 1:SADS_comps)], fullfile(default_path, sprintf("std_SADS_comps_%d.dat", SADS_comps)));
    % writematrix([q DADS(:, 1:SADS_comps)], fullfile(default_path, sprintf("DADS_comps_%d.dat", SADS_comps)));
    % writematrix([q std_DADS(:, 1:SADS_comps)], fullfile(default_path, sprintf("std_DADS_comps_%d.dat", SADS_comps)));
    % 
    % writematrix([profile_SAC(:, 1:SADS_comps)], fullfile(default_path, sprintf("profile_comps_%d.dat", SADS_comps)));
    % writematrix([std_profile_SAC(:, 1:SADS_comps)], fullfile(default_path, sprintf("std_profile_comps_%d.dat", SADS_comps)));

end