clc; clearvars; 
%%
%==================================
% 1. Sample data 
%==================================
default_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuI2_30mM_0002";
% default_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuBr2_150mM_0001";
file_path = "integrate1d only/mask_v1.h5/azim.dat";
% file_path = "integrate1d only/mask_v2.h5/azim.dat";

runs = {
    'scan0015';
    'scan0016';
    'scan0017';
    'scan0018';
    'scan0021';
    'scan0022';
    'scan0023';
    'scan0024';
    'scan0029';
    'scan0030';
    'scan0031';
    'scan0032';
    };

% runs = {
%     'scan0077';
%     'scan0078';
%     'scan0079';
%     'scan0080';
%     'scan0081';
%     };


tds_all = [-3e-9, 0, ...
          100e-12, 178e-12, 316e-12, 562e-12, ...
          1e-9, 1.78e-9, 3.16e-9, 5.62e-9, ...
          10e-9, 17.8e-9, 31.6e-9, 56.2e-9, ...
          100e-9, 178e-9, 316e-9, 562e-9, ...
          1e-6];

% verbose = true;
verbose = false;

save = true;
% save = false;

[mergedData, mergedStd, ~, ~, q] = process_trxl_runs(default_path, runs, file_path, tds_all, [4, 7], true);
if save
    writematrix([q, mergedData], fullfile(default_path, "merged_dat.dat"));
    writematrix([q, mergedStd], fullfile(default_path, "merged_std.dat"));
end

%% SVD on merged data (before PEPCed)
rangeSVD = [1,7];  % SVD 분석용 q 범위
qSVDBool = rangeSVD(1)<q & q<rangeSVD(2);
[U, S, V] = svd(mergedData(qSVDBool, :), 'econ');
if verbose
    HKifuncs.inspect_SVD_v2(U, S, V, tds_all(2:end), 7, 9999, "Merged data (before PEPC)", [], q(qSVDBool, :));
end
%% Plotting: Merged Data
% 시간축 준비: 첫 번째(-3ns)는 Reference이므로 제외하고 2번째부터 사용
if verbose
    plot_times = tds_all(2:end); 
    plot_data_merged = mergedData(qSVDBool, :); % (DataPEPCed 아님)
    plot_q = q(qSVDBool);
    
    % 함수 호출 (q는 SVD 범위 1~7 사용)
    make_qds_plot(plot_q, plot_data_merged, plot_times, 'Merged data (before PEPC)');
    HKifuncs.draw_contour(plot_data_merged, tds_all(2:end), plot_q, 100, "Merged data (before PEPC)", [], [-3e-4 3e-4])
end
%%
%==================================
% 2. Heating data 
%==================================

sol_default_path = "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\heating_MeCN_0001";
% sol_default_path = "C:\Users\IBS\Desktop\TRXL_DG 작업 공간\230425_ESRF_AuBr2\heating_MeCN";
solv_file_path = "integrate1d only/mask_v1.h5/azim.dat";
% solv_file_path = "integrate1d only/mask_v2.h5/azim.dat";
solv_runs = {
    'scan0028';
    'scan0029';
    'scan0030';
    'scan0031';
    'scan0032';
    'scan0033';
    'scan0034';
    'scan0035';
};

td_solv = [-3e-9, 0, ...
          100e-12, 1e-9, 10e-9, 100e-9, ...
          1e-6];
[mergedSolvData, ~, ~, ~, qw] = process_trxl_runs(sol_default_path, solv_runs, solv_file_path, td_solv, [4, 7], true);
if save 
    writematrix([qw, mergedSolvData], fullfile(sol_default_path, "merged_solv_dat.dat"))
end

%%
if verbose
    plot_times = td_solv(2:end); 
    plot_data_solv = mergedSolvData(qSVDBool, :); 
    plot_qw = qw(qSVDBool);
    make_qds_plot(plot_qw, plot_data_solv, plot_times, 'Heating data');
    HKifuncs.draw_contour(plot_data_solv, td_solv(2:end), plot_qw, 100, "Heating data", [], [-5e-4 5e-4])
end

%% SVD on solvent data
rangeSolvSVD = [1,7];  % SVD 분석용 q 범위
solv_qSVDBool = rangeSolvSVD(1)<qw & qw<rangeSolvSVD(2);

[Uw, Sw, Vw] = svd(mergedSolvData(solv_qSVDBool, :), 'econ');
if verbose
    HKifuncs.inspect_SVD_v2(Uw, Sw, Vw, td_solv(2:end), 7, 777, "heating MeCN", [], qw(solv_qSVDBool, :));
end

%%
%==================================
% 3. PEPC 
%==================================
rangePEPC= [1,7];  % PEPC용 q 범위
qPEPC = rangePEPC(1)<q & q<rangePEPC(2);
q_pc  = q(qPEPC);

[~, DataPEPCed] = HKifuncs.pepc(mergedData(qPEPC, :), [Uw(:, 1:4) ones(size(qw(qSVDBool))) 1./qw(qSVDBool)]);

[Up, Sp, Vp] = svd(DataPEPCed, 'econ');
HKifuncs.inspect_SVD_v2(Up, Sp, Vp, tds_all(2:end), 7, 333, "PEPCed data (comps: 4)", [], q_pc);

%% PEPCed data plotting
% 함수 호출 (q 대신 q_pc 사용 필수!)
make_qds_plot(q_pc, DataPEPCed, tds_all(2:end), 'PEPCed data (comps: 4)');
HKifuncs.draw_contour(DataPEPCed, tds_all(2:end), q_pc, 22, "PEPCed data (comps: 4)", [], [-3e-4 3e-4])

%%
if save
    writematrix([q_pc, DataPEPCed], fullfile(default_path, "PEPCed_dat.dat"))
    writematrix([q_pc, mergedStd(qPEPC, :)], fullfile(default_path, "PEPCed_std.dat"))
end

%%
tds_merge = tds_all(2:end) * 1e12;
comps_p = [1 2 3];
[std_p, conc_p] = HKifuncs.fit_V_std(comps_p, DataPEPCed, mergedStd(qPEPC, :), Up, Sp, Vp, 10, tds_merge);

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


%% SADS 뽑아내기
time_constants = xl([10 14 18 19]);
% time_constants = xl([8 11 14 15]);
% time_constants = xl([8 11 12]);


[SAC, std_SAC, theory_profile, profile_SAC, std_profile_SAC] = HKifuncs.KCA(DataPEPCed, mergedStd(qPEPC, :), q(qPEPC), tds_merge, -1, time_constants, tds_merge, 111, [1, 1e6]);
% [DADS, std_DADS, theory_profile_d, profile_SAC_d, std_profile_SAC_d] = HKifuncs.KCA_DADS(data_merge_all, std_merge_all, q, tds_merge, -1, time_constants, tds_merge, 200, [1 1e6]);
%% SADS comp 고른 뒤 저장
SADS_comps = 3;
% DADS_comps = 3;

if save
    writematrix([q(qPEPC) SAC(:, 1:SADS_comps)], fullfile(default_path, sprintf("SADS_comps_%d.dat", SADS_comps)));
    writematrix([q(qPEPC) std_SAC(:, 1:SADS_comps)], fullfile(default_path, sprintf("./std_SADS_comps_%d.dat", SADS_comps)));
    % writematrix([q(qPEPC) DADS(:, 1:DADS_comps)], fullfile(default_path, sprintf("DADS_comps_%d.dat", SADS_comps)));
    % writematrix([q(qPEPC) std_DADS(:, 1:DADS_comps)], fullfile(default_path, sprintf("./std_DADS_comps_%d.dat", SADS_comps)));
end

fprintf('%.6e ', xl(ind));
fprintf("\n");
fprintf('%.6e ', errsl(ind));
%%
%==================================
% 4. 함수들
%==================================
function [mergedData, mergedStd, data4runs, std4runs, q] = process_trxl_runs(default_path, runs, filepath, tds, norm_range, writeFiles)
    if isempty(filepath)
        filepath = 'azim.dat';
    end
    if nargin < 5 || isempty(norm_range)
        norm_range = [4.0 7.0];
    end
    if nargin < 6 || isempty(writeFiles)
        writeFiles = false;
    end
    
    NUM_TIME_DELAY = numel(tds);
    
    data4runs = cell(1, numel(runs));
    std4runs  = cell(1, numel(runs));

    tot_shots = 0; % 전체 샷 수를 카운트하기 위한 변수    q = [];
    for i = 1:numel(runs)

        data_path = fullfile(default_path, runs{i}, filepath);

        dat_all   = readmatrix(data_path);
    
        q        = dat_all(:, 1);
        data     = dat_all(:, 2:end);
        all_shots = size(data, 2);
        
        % 전체 샷 수 누적 (나중에 sqrt(N)으로 나눌 때 사용)
        tot_shots = tot_shots + all_shots;
        [~, data] = DHanfuncs.normalize_AUC(q, data, norm_range);
    
        NUM_SHOTS  = all_shots / NUM_TIME_DELAY;
        if mod(all_shots, NUM_TIME_DELAY) ~= 0
            error('Error: 파일 %s의 데이터 열 개수가 시간 지연 개수의 정수배가 아닙니다.', runs{i});
        end

        split_size = repmat(NUM_TIME_DELAY, 1, NUM_SHOTS);
    
        chunk_data = mat2cell(data, size(data, 1), split_size);
    
        chunk_diff_data = cell(1, NUM_SHOTS);
        for shot = 1:NUM_SHOTS
            data4shot           = chunk_data{shot};
            chunk_diff_data{shot} = data4shot(:, 2:end) - data4shot(:, 1);
        end
    
        diff_data = cat(3, chunk_diff_data{:});
        avgData   = mean(diff_data, 3);
        stdData   = std(diff_data, 0, 3);
    
        data4runs{i} = avgData;
        std4runs{i}  = stdData;

        if writeFiles
            folderName = fullfile(default_path, runs{i}, 'DiffAve');
            if ~isfolder(folderName)
                mkdir(folderName);
            end
    
            for td = 1:NUM_TIME_DELAY-1  % negative td는 뺀다.
                filename = sprintf('diff_av_%.6e.txt', tds(td+1));  % pos td부터
                outpath  = fullfile(folderName, filename);
                writematrix([q avgData(:, td) stdData(:, td) q.*avgData(:, td)], ...
                            outpath, 'Delimiter', '\t');
            end
        end
    end
    
    % 1. Cell Array를 3차원 행렬로 변환 [q, delay, run]
    all_means = cat(3, data4runs{:}); % data_all에 해당
    all_stds  = cat(3, std4runs{:});  % std_all에 해당
    
    mergedData = mean(all_means, 3);
    mergedStd  = mean(all_stds, 3);
    mergedStd  = mergedStd ./ sqrt(tot_shots);


    % % % % 2. 가중치(샷 수) 준비: 계산을 위해 3차원으로 reshape [1, 1, num_runs]
    % % % W = reshape(shots_per_run, 1, 1, []);  % 1, 1, 12 size
    % % % sum_W = sum(W, 3); % 총 샷 수 (분모용)
    % % % 
    % % % % 3. Merged Data (가중 평균)
    % % % % 식: sum(N_i * Mean_i) / sum(N_i)
    % % % mergedData = sum(W .* all_means, 3) ./ sum_W;
    % % % 
    % % % % 4. Merged Std (복합 오차 공식)
    % % % % 식: sqrt( [sum(N^2 * std^2) + sum(N * mean^2)] / (sum_N)^2 - (merged_mean^2 / sum_N) )
    % % % 
    % % % term1 = sum( (W.^2) .* (all_stds.^2), 3 );      % 내부 분산 항
    % % % term2 = sum( W .* (all_means.^2), 3 );          % Run 간 편차 항
    % % % 
    % % % temp_sq_sum = term1 + term2;
    % % % 
    % % % mergedStd = sqrt( (temp_sq_sum ./ (sum_W.^2)) - ((mergedData.^2) ./ sum_W) );
    
end




function make_qds_plot(q_vec, data_mat, time_vec, fig_title)
% MAKE_WATERFALL_PLOT Waterfall 형태의 시분해 산란 그래프를 그리는 함수
%
%   q_vec:      q 축 벡터 (N x 1)
%   data_mat:   데이터 행렬 (N x TimePoints) - MergedData 또는 DataPEPCed
%   time_vec:   시간 지연 벡터 (1 x TimePoints) - 데이터 열 개수와 같아야 함
%   fig_title:  그래프 제목 (문자열)

    % 1. 예외 처리: q와 데이터 길이 확인
    if length(q_vec) ~= size(data_mat, 1)
        error('q 벡터의 길이와 데이터의 행(row) 개수가 맞지 않습니다.');
    end
    if length(time_vec) ~= size(data_mat, 2)
        error('시간 벡터의 길이와 데이터의 열(column) 개수가 맞지 않습니다.');
    end

    % 2. 시간 오름차순 정렬 (Sorting)
    [sorted_times, sort_idx] = sort(time_vec, 'ascend');
    sorted_data = data_mat(:, sort_idx);

    % 3. q-weighting (q * dS)
    plot_data = sorted_data .* q_vec;

    % 4. 오프셋 및 색상 설정
    num_curves = length(sorted_times);
    max_amp = max(abs(plot_data(:)));
    offset_step = max_amp * 0.6; % 간격 조절

    curve_colors = cool(num_curves); % 색상 (Cyan -> Magenta)

    % 5. 그래프 그리기
    figure('Color', 'w', 'Position', [100, 100, 900, 800], 'Name', fig_title);
    hold on;
    
    x_text_pos = max(q_vec); % 라벨 달 위치 (오른쪽 끝)

    for i = 1:num_curves
        current_y = plot_data(:, i);
        current_t = sorted_times(i);
        
        % 오프셋 계산 (위에서 아래로)
        current_offset = -(i-1) * offset_step;
        
        % (A) Zero Line (기준선)
        yline(current_offset, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
        % (B) 데이터 플롯
        plot(q_vec, current_y + current_offset, 'LineWidth', 2, 'Color', curve_colors(i, :));
        
        % (C) 시간 라벨 (기준선 높이, 중앙 정렬)
        label_str = local_format_time(current_t);
        text(x_text_pos, current_offset, ['  ' label_str], ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
    end

    % 6. 스타일 꾸미기
    xlabel('q', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('q \Delta S(q)', 'FontSize', 14, 'FontWeight', 'bold');
    title(fig_title, 'FontSize', 20, 'FontWeight', 'bold');
    
    xlim([min(q_vec), max(q_vec)]);
    % set(gca, 'YTick', [], 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');
    hold off;
end

% (내부 함수) 시간 단위 변환
function str = local_format_time(t)
    if t == 0
        str = '0 ps';
    elseif abs(t) < 1e-9
        str = sprintf('%.0f ps', t * 1e12);
    elseif abs(t) < 1e-6
        val_ns = t * 1e9;
        if abs(val_ns - round(val_ns)) < 1e-3
            str = sprintf('%.0f ns', val_ns);
        else
            str = sprintf('%.2g ns', val_ns);
        end
    else
        str = sprintf('%.2g \\mus', t * 1e6);
    end
end
