clc; clearvars; 
%%
%==================================
% 1. Sample data 
%==================================
default_path =  "\\172.30.150.180\homes\sdlab\230425_ESRF_AuBr2\SCRIPTS\inHouseProcess\resultsCD\AuI2_30mM_0002";
file_path = "integrate1d only/mask_v1.h5/azim.dat";
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

tds_all = [-3e-9, 0, ...
          100e-12, 178e-12, 316e-12, 562e-12, ...
          1e-9, 1.78e-9, 3.16e-9, 5.62e-9, ...
          10e-9, 17.8e-9, 31.6e-9, 56.2e-9, ...
          100e-9, 178e-9, 316e-9, 562e-9, ...
          1e-6];
save = true;
% save = false;

[mergedData, mergedStd, ~, ~, q] = process_trxl_runs_weighted(default_path, runs, file_path, tds_all, [4, 7], true);
if save
    writematrix([q, mergedData], fullfile(default_path, "merged_dat_v2.dat"));
    writematrix([q, mergedStd], fullfile(default_path, "merged_std_v2.dat"));
end
%% SVD on merged data (before PEPCed)
rangeSVD = [1,7];  % SVD 분석용 q 범위
qSVDBool = rangeSVD(1)<q & q<rangeSVD(2);

[U, S, V] = svd(mergedData(qSVDBool, :), 'econ');
HKifuncs.inspect_SVD_v2(U, S, V, tds_all(2:end), 7, 9999, "Merged data (before PEPC)", [1, 5], q(qSVDBool, :));

%% Plotting: Merged Data
% 시간축 준비: 첫 번째(-3ns)는 Reference이므로 제외하고 2번째부터 사용
plot_times = tds_all(2:end); 
plot_data_merged = mergedData(qSVDBool, :); % (DataPEPCed 아님)
plot_q = q(qSVDBool);

% 함수 호출 (q는 SVD 범위 1~7 사용)
make_qds_plot(plot_q, plot_data_merged, plot_times, 'Merged data (before PEPC)');
HKifuncs.draw_contour(plot_data_merged, tds_all(2:end), plot_q, 100, "Merged data (before PEPC)", [], [-2e-3 2e-3])

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
% solv_runs = {
%     'run0028';
%     'run0029';
%     'run0030';
%     'run0031';
%     'run0032';
%     'run0033';
%     'run0034';
%     'run0035';
% };
td_solv = [-3e-9, 0, ...
          100e-12, 1e-9, 10e-9, 100e-9, ...
          1e-6];
[mergedSolvData, ~, ~, ~, qw] = process_trxl_runs_weighted(sol_default_path, solv_runs, solv_file_path, td_solv, [4, 7], true);
if save 
    writematrix([qw, mergedSolvData], fullfile(sol_default_path, "merged_solv_dat_v2.dat"))
end
%% SVD on solvent data
rangeSolvSVD = [1,7];  % SVD 분석용 q 범위
solv_qSVDBool = rangeSolvSVD(1)<qw & qw<rangeSolvSVD(2);

[Uw, Sw, Vw] = svd(mergedSolvData(solv_qSVDBool, :), 'econ');
HKifuncs.inspect_SVD_v2(Uw, Sw, Vw, td_solv(2:end), 7, 777, "heating MeCN", [], qw(solv_qSVDBool, :));

%%
%==================================
% 3. PEPC 
%==================================
rangePEPC= [1,7];  % PEPC용 q 범위
qPEPC = rangePEPC(1)<q & q<rangePEPC(2);
q_pc  = q(qPEPC);

[~, DataPEPCed] = HKifuncs.pepc(mergedData(qPEPC, :), Uw(:, 1:3));

[Up, Sp, Vp] = svd(DataPEPCed, 'econ');
HKifuncs.inspect_SVD_v2(Up, Sp, Vp, tds_all(2:end), 7, 333, "PEPCed data (comps: 3)", [], q_pc);

%% PEPCed data plotting
% 함수 호출 (q 대신 q_pc 사용 필수!)
make_qds_plot(q_pc, DataPEPCed, tds_all(2:end), 'PEPCed data (comps: 3)');
HKifuncs.draw_contour(DataPEPCed, tds_all(2:end), q_pc, 22, "PEPCed data (comps: 3)", [], [-5e-4 5e-4])

%%
if save
    writematrix([q_pc, DataPEPCed], fullfile(default_path, "PEPCed_dat_v2.dat"))
    writematrix([q_pc, mergedStd(qPEPC, :)], fullfile(default_path, "PEPCed_std_v2.dat"))
end
%%

function [mergedData, mergedStd, data4runs, std4runs, q, tds_merge] = process_trxl_runs_weighted(default_path, runs, file_path, tds_merge, norm_range, writeDiffAve)
% PROCESS_TRXL_RUNS_WEIGHTED 
%   샷(Shot) 수를 가중치(w)로 사용하여 여러 Run을 병합하는 함수
%   
%   입력:
%       tds_merge: 최종적으로 병합할 목표 시간 축 (Target Time Grid)
    if isempty(file_path)
        file_path = 'azim.dat';
    end
    if nargin < 5 || isempty(norm_range)
        norm_range = [4.0 7.0];
    end
    if nargin < 6 || isempty(writeDiffAve)
        writeDiffAve = false;
    end
    
    % 1. 데이터를 모을 컨테이너 초기화
    % (나중에 binning 알고리즘에 넣기 위해 모든 Run의 데이터를 옆으로 쭉 붙입니다)
    grand_data = [];  % data_all 역할
    grand_std  = [];  % std_all 역할
    grand_w    = [];  % imgs_all (가중치) 역할
    grand_td   = [];  % td_all 역할
    
    data4runs = cell(1, numel(runs));
    std4runs  = cell(1, numel(runs));
    q = [];
    
    % =========================================================
    % 2. 각 Run 파일 로드 및 1차 가공 (샷 단위 평균/표준편차)
    % =========================================================
    for i = 1:numel(runs)
        data_path = fullfile(default_path, runs{i}, file_path);
        dat_all   = readmatrix(data_path);
    
        q_current = dat_all(:, 1);
        tds_raw   = dat_all(1, 2:end);     % 이 파일에 기록된 전체 Time points (반복 포함)
        data      = dat_all(:, 2:end);
        
        % q 벡터 저장 (첫 번째 Run 기준)
        if isempty(q)
            q = q_current;
        elseif ~isequal(q, q_current)
            % (Q축이 다르면 보간이 필요하지만 여기선 에러 처리)
            error('Run마다 q축이 다릅니다. 확인이 필요합니다.');
        end
    
        % 정규화 (Normalization)
        normfactor = AUC_Jkim(data, q, norm_range);
        data       = data ./ normfactor;
    
        % 샷(Shot) 분리 및 통계 계산
        % tds_merge(입력받은 시간축)의 길이를 기준으로 샷 개수 파악
        NUM_TIME_DELAY = numel(tds_merge); 
        NUM_SHOTS      = numel(tds_raw) / NUM_TIME_DELAY;
        
        % 데이터 무결성 체크
        if mod(NUM_SHOTS, 1) ~= 0
            error('%s: 데이터 열 개수가 시간 축 개수의 정수배가 아닙니다.', runs{i});
        end
        
        split_size = repmat(NUM_TIME_DELAY, 1, NUM_SHOTS);
        chunk_data = mat2cell(data, size(data, 1), split_size);
    
        chunk_diff_data = cell(1, NUM_SHOTS);
        for shot = 1:NUM_SHOTS
            data4shot = chunk_data{shot};
            % (Laser On - Laser Off) 계산: 첫 번째 time point를 Reference(Off)로 가정
            chunk_diff_data{shot} = data4shot(:, 2:end) - data4shot(:, 1);
        end
        
        % [3차원 배열]: (q x TimePoints-1 x Shots)
        diff_data_3d = cat(3, chunk_diff_data{:});
        
        % Run 별 평균 및 표준편차 계산
        avgData_run = mean(diff_data_3d, 3);
        stdData_run = std(diff_data_3d, 0, 3);
        
        data4runs{i} = avgData_run;
        std4runs{i}  = stdData_run;
        
        % ---------------------------------------------------------
        % [Grand Accumulation] 병합을 위해 데이터 누적
        % ---------------------------------------------------------
        % 이 Run에서 사용된 시간축 (Reference인 첫번째 타임 제외)
        % 입력받은 tds_merge의 2번째부터 끝까지가 실제 데이터 시간이라고 가정
        current_tds = tds_merge(2:end); 
        
        % 가중치(w): 이 Run의 샷 수 (모든 시간 포인트에 동일하게 적용)
        current_w   = repmat(NUM_SHOTS, 1, numel(current_tds));
        
        grand_data = [grand_data, avgData_run]; % 옆으로 붙임
        grand_std  = [grand_std, stdData_run];
        grand_w    = [grand_w, current_w];
        grand_td   = [grand_td, current_tds];
        
        % (개별 파일 저장 옵션은 유지)
        if writeDiffAve
            folderName = fullfile(default_path, runs{i}, 'DiffAve');
            if ~isfolder(folderName)
                mkdir(folderName);
            end
            for t_idx = 1:length(current_tds)
                filename = sprintf('diff_av_%.6e.txt', current_tds(t_idx));
                outpath  = fullfile(folderName, filename);
                writematrix([q avgData_run(:, t_idx) stdData_run(:, t_idx) q.*avgData_run(:, t_idx)], ...
                            outpath, 'Delimiter', '\t');
            end
        end
    end
    
    % =========================================================
    % 3. 가중 평균 병합 (Weighted Merging) - 제공해주신 알고리즘 적용
    % =========================================================
    
    % 변수 이름 매핑 (제공해주신 코드와 맞춤)
    data_all = grand_data;
    std_all  = grand_std;
    imgs_all = grand_w;     % 가중치 = 샷 수
    td_all   = grand_td;    
    
    % 병합할 타겟 시간축 (Reference 0 제외)
    target_tds = tds_merge(2:end); 
    
    % --- [User Snippet Start] ---
    
    num_time_bins = numel(target_tds);                   
    num_q         = size(data_all, 1);
    
    mergedData    = zeros(num_q, num_time_bins);
    mergedStd     = zeros(num_q, num_time_bins);
    
    % 최근접 bin 매핑 (각 td_all → 최근접 target_tds 인덱스)
    [~, neareas_bin_index] = min(abs(td_all(:) - target_tds(:).'), [], 2);
    
    for bin_idx = 1:num_time_bins
        target_cols = (neareas_bin_index == bin_idx);
        
        if ~any(target_cols)
            fprintf('Warning: Time bin %g has no data.\n', target_tds(bin_idx));
            continue; 
        end
        
        w    = imgs_all(target_cols);  % 가중치(스캔 수)
        Wsum = sum(w);
        
        % 가중 평균
        % data_all(:, target_cols) : [q x k]
        % w : [1 x k]
        % (w * data_all')' => (1 x q)' => (q x 1)
        mergedData(:, bin_idx) = (w * data_all(:, target_cols)')' / Wsum;
        
        % 오차 전파: sqrt( sum(w^2 * sigma^2) ) / sum(w)
        s2 = std_all(:, target_cols).^2;
        % w.^2 .* s2 : w제곱을 각 q행마다 곱함 (Implicit expansion)
        mergedStd(:, bin_idx) = sqrt( sum( (w.^2) .* s2, 2 ) ) / Wsum;
    end
end

function AUCs = AUC_Jkim(sqMatrix,q,ROI)
    qTF = ROI(1) <= q & q <= ROI(2);
    AUCs = trapz(q(qTF),sqMatrix(qTF,:));
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
