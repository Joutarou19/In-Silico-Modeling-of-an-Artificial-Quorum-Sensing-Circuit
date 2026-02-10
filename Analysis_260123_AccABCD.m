%% Hill式フィッティング：Control vs AccABCD プロフェッショナル比較解析
clear; clc; close all;

%% 1. 設定
filename = 'C:\Users\jouta\Downloads\Data_AccABCD挿入株.xlsx';
sheetname = 'GFP-0hGFP%OD';

% 解析対象の定義（名前, 開始行, 終了行, 配色RGB）
strain_info = {
    'Control', 100, 107, [0 0.447 0.741];  % 鮮やかな青
    'AccABCD', 113, 120, [0.85 0.325 0.098] % 鮮やかなオレンジ
};

%% 2. 最適化エンジンの設定
hill_model = @(p, x) (p(1) .* x.^p(3)) ./ (p(2) + x.^p(3));
lsq_options = optimoptions('lsqcurvefit', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
ms = MultiStart('Display', 'off');
results = struct();

% グラフの土台作成
figure('Name', 'Sensitivity Comparison', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
hold on;

%% 3. ループ処理による各株の解析
plot_handles = []; % 凡例用ハンドル

for i = 1:size(strain_info, 1)
    name = strain_info{i, 1};
    r_start = strain_info{i, 2};
    r_end = strain_info{i, 3};
    color = strain_info{i, 4};
    
    % --- データ読み込み ---
    iptg_range = sprintf('B%d:B%d', r_start, r_end);
    data_range = sprintf('C%d:EJ%d', r_start, r_end);
    iptg_concs = readmatrix(filename, 'Sheet', sheetname, 'Range', iptg_range);
    raw_data = readmatrix(filename, 'Sheet', sheetname, 'Range', data_range);
    max_response = max(raw_data, [], 2, 'omitnan');
    
    % --- 最適化実行 ---
    v0 = max(max_response);
    initial_p = [v0, (mean(iptg_concs))^2, 2];
    problem = createOptimProblem('lsqcurvefit', 'objective', hill_model, 'x0', initial_p, ...
        'xdata', iptg_concs, 'ydata', max_response, 'lb', [0, 0, 0.1], 'ub', [v0*5, 10, 20], 'options', lsq_options);
    [p_fit, ~] = run(ms, problem, 30);
    
    % 結果の算出
    ec50 = p_fit(2)^(1/p_fit(3));
    results.(name).p = p_fit;
    results.(name).ec50 = ec50;
    
    % --- プロット ---
    % フィッティング曲線
    x_smooth = logspace(-3, 1.5, 500); % 範囲を少し広めに設定
    y_fit = hill_model(p_fit, x_smooth);
    h_fit = semilogx(x_smooth, y_fit, '-', 'Color', color, 'LineWidth', 3, ...
        'DisplayName', sprintf('%s (EC_{50}: %.3f mM)', name, ec50));
    
    % 実験データ（白抜き＋色枠のマーカーで視認性を向上）
    h_data = semilogx(iptg_concs, max_response, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off'); % 凡例は曲線に任せる
    
    % EC50位置の垂直線（ドット線）
    xline(ec50, ':', 'Color', color, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    plot_handles = [plot_handles, h_fit];
end

%% 4. グラフの装飾（右図のスタイルに合わせる）
grid on; 
ax = gca;
ax.GridLineStyle = '-';
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.XScale = 'log';
ax.XLim = [1e-3, 1e1]; % 表示範囲の調整
box on;

% フォント設定
set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
xlabel('IPTG Concentration (mM)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Max GFP/OD_{600}', 'FontSize', 16, 'FontWeight', 'bold');
title('Sensitivity Comparison (Hill Equation Fit)', 'FontSize', 18);

% 凡例の改善（左上に配置）
lgd = legend(plot_handles, 'Location', 'northwest', 'FontSize', 13);
title(lgd, 'Strain & Sensitivity');
lgd.EdgeColor = [0.5 0.5 0.5];

% 右図のような詳細凡例（右下）が必要な場合は追加
% annotation('textbox', [0.65 0.15 0.25 0.1], 'String', {'● Control (Data)', '● AccABCD (Data)'}, ...
%     'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', [0.8 0.8 0.8]);

%% 5. 結果レポート出力
fprintf(['\n', repmat('=',1,50), '\n']);
fprintf('  最終解析レポート\n');
fprintf([repmat('-',1,50), '\n']);
f = fieldnames(results);
for i = 1:length(f)
    n = f{i};
    fprintf('【%s】 EC50: %.4f mM | Hill n: %.2f | Vmax: %.1f\n', ...
        n, results.(n).ec50, results.(n).p(3), results.(n).p(1));
end
fprintf([repmat('=',1,50), '\n']);