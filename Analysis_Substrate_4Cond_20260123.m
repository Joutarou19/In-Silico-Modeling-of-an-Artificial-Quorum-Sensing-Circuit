%% 初期設定
clear; clc; close all;

% ファイルパスの設定
filename = 'C:\Users\jouta\OneDrive\デスクトップ\研究\トレーニング\Data_基質添加実験.xlsx'; % 必要に応じてフルパスに変更してください
sheetname = 'GFP-0hGFP%OD';

%% 1. データの読み込み
% ファイル構造に基づき、各条件のデータ範囲を指定します
try
    % --- Control (N) ---
    % IPTG: C103:C110, Data: E103:DT110
    ctrl_iptg = readmatrix(filename, 'Sheet', sheetname, 'Range', 'C103:C110');
    ctrl_data = readmatrix(filename, 'Sheet', sheetname, 'Range', 'E103:DT110');
    
    % --- Malonate ---
    % IPTG: C111:C118, Data: E111:DT118
    malo_iptg = readmatrix(filename, 'Sheet', sheetname, 'Range', 'C111:C118');
    malo_data = readmatrix(filename, 'Sheet', sheetname, 'Range', 'E111:DT118');
    
    % --- Pyruvate ---
    % IPTG: C119:C126, Data: E119:DT126
    pyru_iptg = readmatrix(filename, 'Sheet', sheetname, 'Range', 'C119:C126');
    pyru_data = readmatrix(filename, 'Sheet', sheetname, 'Range', 'E119:DT126');
    
    % --- Methionine (MetK) ---
    % IPTG: C127:C134, Data: E127:DT134
    met_iptg  = readmatrix(filename, 'Sheet', sheetname, 'Range', 'C127:C134');
    met_data  = readmatrix(filename, 'Sheet', sheetname, 'Range', 'E127:DT134');
    
catch ME
    error('データの読み込みに失敗しました。\nエラー: %s', ME.message);
end

%% 2. データ前処理（最大値の取得）
% 各条件について、定常状態（最大応答）を計算
ctrl_resp = max(ctrl_data, [], 2, 'omitnan');
malo_resp = max(malo_data, [], 2, 'omitnan');
pyru_resp = max(pyru_data, [], 2, 'omitnan');
met_resp  = max(met_data,  [], 2, 'omitnan');

% データを構造体にまとめる
conditions(1).name = 'Control';    conditions(1).iptg = ctrl_iptg; conditions(1).resp = ctrl_resp; conditions(1).color = 'b';
conditions(2).name = 'Malonate';   conditions(2).iptg = malo_iptg; conditions(2).resp = malo_resp; conditions(2).color = [0 0.5 0]; % Dark Green
conditions(3).name = 'Pyruvate';   conditions(3).iptg = pyru_iptg; conditions(3).resp = pyru_resp; conditions(3).color = 'm'; % Magenta
conditions(4).name = 'Methionine'; conditions(4).iptg = met_iptg;  conditions(4).resp = met_resp;  conditions(4).color = 'r';

%% 3. Hill方程式によるフィッティング (厳密なEC50算出)
% モデル: Y = Bottom + (Top - Bottom) * [S]^n / (Kd + [S]^n)
% p = [Top, Kd, n, Bottom]
hill_model = @(p, x) p(4) + (p(1) - p(4)) .* (x.^p(3)) ./ (p(2) + x.^p(3));

% 目的関数
obj_func = @(p, x, y) sum((y - hill_model(p, x)).^2, 'omitnan');

% 最適化オプション
opts = optimset('Display', 'off', 'MaxFunEvals', 20000, 'MaxIter', 20000);

fprintf('--- 解析結果 (EC50 = Kd^(1/n)) ---\n');

for i = 1:length(conditions)
    x = conditions(i).iptg;
    y = conditions(i).resp;
    
    % 初期値推定
    top_0 = max(y);
    bot_0 = min(y);
    n_0 = 2;
    kd_0 = (0.1)^n_0;
    
    p0 = [top_0, kd_0, n_0, bot_0];
    
    % フィッティング実行
    p_fit = fminsearch(@(p) obj_func(p, x, y), p0, opts);
    
    % パラメータ抽出とEC50算出
    Top = p_fit(1);
    Kd = p_fit(2);
    n = p_fit(3);
    Bottom = p_fit(4);
    EC50 = Kd^(1/n); 
    
    % 結果を保存
    conditions(i).p_fit = p_fit;
    conditions(i).EC50 = EC50;
    conditions(i).n = n;
    
    fprintf('%-12s: EC50 = %.4f mM (n = %.2f, Max = %.0f)\n', ...
        conditions(i).name, EC50, n, Top);
end
fprintf('----------------------------------\n');

%% 4. グラフ描画 (対数スケールのみ)
figure('Name', 'Sensitivity Comparison (Log Only)', 'Position', [100, 100, 1000, 600], 'Color', 'w');
hold on;

% プロット用X軸作成
x_smooth_log = logspace(log10(0.04), log10(0.9), 300);

for i = 1:length(conditions)
    % 実測値プロット
    scatter(conditions(i).iptg, conditions(i).resp, 80, conditions(i).color, ...
        'filled', 'MarkerEdgeColor', 'k', 'DisplayName', conditions(i).name);
    
    % 近似曲線プロット
    y_fit = hill_model(conditions(i).p_fit, x_smooth_log);
    plot(x_smooth_log, y_fit, '-', 'Color', conditions(i).color, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    % EC50ライン
    xline(conditions(i).EC50, ':', 'Color', conditions(i).color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % EC50値のテキスト表示 (左上にまとめて表示)
    text(0.05, 0.95 - (i-1)*0.06, sprintf('%s EC_{50}: %.3f mM', conditions(i).name, conditions(i).EC50), ...
        'Units', 'normalized', 'Color', conditions(i).color, 'FontWeight', 'bold', 'FontSize', 12, ...
        'BackgroundColor', 'w', 'EdgeColor', 'none');
end

% 軸設定
set(gca, 'XScale', 'log');
xlim([0.04, 0.9]); % X軸の範囲はデータの分布に合わせて適宜調整してください
title('Sensitivity Comparison (Log Scale)', 'FontSize', 16);
xlabel('IPTG Concentration (mM)', 'FontSize', 14);
ylabel('Max GFP/OD_{600}', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 12);

% グリッドと枠線

grid on; grid minor; box on;
hold off;

plotedit on; 
inspect(gcf);  