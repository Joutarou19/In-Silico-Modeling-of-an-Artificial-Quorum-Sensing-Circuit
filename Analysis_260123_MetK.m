%% Hill式フィッティング：論文品質 比較解析（アドオン活用版）
clear; clc; close all;

%% 1. 設定
filename = 'C:\Users\jouta\OneDrive\デスクトップ\研究\Original_Data\Data_MetK挿入株.xlsx';
sheetname = '(GFP-0hGFP)%OD';

% 解析対象の定義（名前, 行範囲, [R G B]配色）
strains = {
    'Control', [100, 111], [0 0.447 0.741];   % 鮮やかな青
    'MetK',    [113, 124], [0.85 0.325 0.098]  % 鮮やかなオレンジ/赤
};

%% 2. 最適化エンジンの設定
hill_model = @(p, x) (p(1) .* x.^p(3)) ./ (p(2) + x.^p(3));
lsq_opts = optimoptions('lsqcurvefit', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
ms = MultiStart('Display', 'iter');
results = struct();

% グラフの土台作成（少し横長に設定）
figure('Name', 'Professional Sensitivity Comparison', 'Color', 'w', 'Position', [100, 100, 950, 600]);
hold on;

%% 3. 解析ループ
line_handles = []; % 凡例用

for i = 1:size(strains, 1)
    name  = strains{i, 1};
    rows  = strains{i, 2};
    color = strains{i, 3};
    
    fprintf('\n>>> %s 株の解析を実行中...\n', name);
    
    % --- データ読み込み ---
    iptg_range = sprintf('B%d:B%d', rows(1), rows(2));
    data_range = sprintf('C%d:EJ%d', rows(1), rows(2));
    iptg_raw = readmatrix(filename, 'Sheet', sheetname, 'Range', iptg_range);
    data_raw = readmatrix(filename, 'Sheet', sheetname, 'Range', data_range);
    
    valid = ~isnan(iptg_raw);
    x_data = iptg_raw(valid);
    y_data = max(data_raw(valid, :), [], 2, 'omitnan');
    
    % --- 最適化実行 ---
    v0 = max(y_data);
    problem = createOptimProblem('lsqcurvefit', 'objective', hill_model, ...
        'x0', [v0, (mean(x_data))^2, 2], 'xdata', x_data, 'ydata', y_data, ...
        'lb', [0, 0, 0.1], 'ub', [v0*5, 10, 20], 'options', lsq_opts);
    [p_fit, ~] = run(ms, problem, 30);
    
    % 結果算出
    ec50 = p_fit(2)^(1/p_fit(3));
    results.(name).p = p_fit;
    results.(name).ec50 = ec50;
    
    % --- プロット ---
    % 1. フィッティング曲線
    x_smooth = logspace(log10(1e-3), log10(max(x_data)*2), 500);
    y_fit = hill_model(p_fit, x_smooth);
    h_fit = semilogx(x_smooth, y_fit, '-', 'Color', color, 'LineWidth', 3, ...
        'DisplayName', sprintf('%s EC_{50}: %.3f mM', name, ec50));
    
    % 2. 実験データ（白枠を付けて視認性をアップ）
    semilogx(x_data, y_data, 'o', 'MarkerSize', 9, 'MarkerFaceColor', color, ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    % 3. EC50垂直線（ドット線）
    xline(ec50, ':', 'Color', color, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    line_handles = [line_handles, h_fit];
end

%% 4. グラフ装飾（画像のデザインを再現）
grid on; box on;
ax = gca;
ax.XScale = 'log';
ax.FontSize = 14;
ax.LineWidth = 1.2;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.1;
ax.XMinorGrid = 'on';

% 軸ラベルとタイトル
xlabel('IPTG Concentration (mM)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Max GFP/OD_{600}', 'FontSize', 16, 'FontWeight', 'bold');
title('Sensitivity Comparison (Hill Equation Fit)', 'FontSize', 18);

% 凡例の改善（左上に配置）
lgd = legend(line_handles, 'Location', 'northwest', 'FontSize', 13);
title(lgd, 'Strain & Sensitivity');
lgd.EdgeColor = [0.6 0.6 0.6]; % 凡例の枠を少し薄く

% 余白の調整
xlim([1e-3, 1e1]);

%% 5. レポート表示
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