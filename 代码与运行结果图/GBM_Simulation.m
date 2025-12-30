%% 应用随机过程大作业：几何布朗运动 (GBM) 模拟与分析
% 包含：轨道生成、参数敏感性分析、二次变差验证、实证参数估计
clear; clc; close all;
rng(2025); % 设置随机种子，保证每次运行结果一致

%% 1. 模型参数设置 (Basic Settings)
mu = 0.15;      % 漂移率 (期望收益率)
sigma = 0.30;   % 波动率 (风险)
S0 = 100;       % 初始价格
T = 1.0;        % 时间长度 (1年)
dt = 0.001;     % 时间步长
N = floor(T/dt);% 总步数
t = (0:dt:T)';  % 时间向量

%% 2. 模拟多条轨道 (Simulation of Paths) - 对应要求(c)
num_paths = 10; % 模拟轨道的数量

% 预分配矩阵
S = zeros(N+1, num_paths);
S(1, :) = S0;

% 生成标准布朗运动增量 dW ~ N(0, dt)
dW = sqrt(dt) * randn(N, num_paths);

% 方法：使用几何布朗运动的解析解 (Exact Solution)
% S(t) = S0 * exp( (mu - 0.5*sigma^2)t + sigma*W(t) )
W = [zeros(1, num_paths); cumsum(dW)]; % 布朗运动路径
time_matrix = repmat(t, 1, num_paths); % 时间矩阵扩充

% 计算股价路径
S = S0 .* exp((mu - 0.5 * sigma^2) .* time_matrix + sigma .* W);

% --- 绘图 1: GBM 多条轨道 ---
figure('Name', 'GBM Simulation', 'Color', 'w');
plot(t, S, 'LineWidth', 1.2);
xlabel('Time (Year)');
ylabel('Price ($)');
title(['几何布朗运动模拟 (Paths = ', num2str(num_paths), ')']);
grid on;
hold on;
plot(t, mean(S, 2), 'k--', 'LineWidth', 2); % 绘制均值路径
legend('Simulated Paths', 'Mean Path', 'Location', 'northwest');

%% 3. 参数敏感性分析 (Parameter Sensitivity) - 对应要求(d)
% 探究不同波动率 sigma 对轨道的影响

sigma_low = 0.1;
sigma_high = 0.6;

% 计算两种不同波动率下的路径 (为了对比，使用相同的随机种子/增量)
S_low  = S0 .* exp((mu - 0.5 * sigma_low^2) .* t + sigma_low .* W(:,1));
S_high = S0 .* exp((mu - 0.5 * sigma_high^2) .* t + sigma_high .* W(:,1));

% --- 绘图 2: 参数对比 ---
figure('Name', 'Sensitivity Analysis', 'Color', 'w');
plot(t, S_low, 'b-', 'LineWidth', 1.5); hold on;
plot(t, S_high, 'r-', 'LineWidth', 1.5);
xlabel('Time (Year)');
ylabel('Price');
title('波动率参数对轨道的影响 (\sigma_{low}=0.1 vs \sigma_{high}=0.6)');
legend(['Low Volatility (\sigma=', num2str(sigma_low), ')'], ...
       ['High Volatility (\sigma=', num2str(sigma_high), ')']);
grid on;

%% 4. 二次变差验证 (Quadratic Variation) - 对应要求(d)
% 理论：对数收益率的二次变差应收敛于 sigma^2 * T
% 二次变差 = sum( (log(S_{t+1}) - log(S_t))^2 )

% 提取第一条轨道进行分析
S_path = S(:, 1);
log_S = log(S_path);
log_returns = diff(log_S); % 对数差分

% 计算累积二次变差
realized_QV = cumsum(log_returns.^2);
theoretical_QV = (sigma^2) * t(2:end); % 理论值 sigma^2 * t

% --- 绘图 3: 二次变差 ---
figure('Name', 'Quadratic Variation', 'Color', 'w');
plot(t(2:end), realized_QV, 'b-', 'LineWidth', 1.5); hold on;
plot(t(2:end), theoretical_QV, 'r--', 'LineWidth', 2);
xlabel('Time');
ylabel('Cumulative Quadratic Variation');
title('二次变差收敛性验证');
legend('Realized QV (Simulated)', 'Theoretical QV (\sigma^2 t)', 'Location', 'northwest');
grid on;
% 结论：如果两条线重合度高，证明了伊藤积分性质

%% 5. 实证分析部分 (Empirical Analysis) - 对应要求(d)(e)
% 说明：实际操作时，请使用 readtable 或 xlsread 读取真实的 CSV 股票数据
% 这里为了代码可直接运行，我们生成一段“伪造”的真实数据

% --- 5.1 模拟读取数据 ---
% 假设这是你从 Yahoo Finance 下载的某只股票一年的收盘价
real_data_len = 252; % 一年交易日
real_prices = 100 + cumsum(randn(real_data_len, 1)); % 随机游走生成伪数据
real_time = linspace(0, 1, real_data_len)';

% --- 5.2 参数估计 ---
% 计算对数收益率
r_real = diff(log(real_prices));
dt_real = 1/252; % 每日的时间步长

% 估计波动率 sigma_hat (标准差 / sqrt(dt))
sigma_hat = std(r_real) / sqrt(dt_real);

% 估计漂移率 mu_hat (均值/dt + 0.5*sigma^2)
mu_hat = mean(r_real) / dt_real + 0.5 * sigma_hat^2;

fprintf('--- 实证分析结果 ---\n');
fprintf('估计的年化波动率 (Sigma): %.4f\n', sigma_hat);
fprintf('估计的年化漂移率 (Mu):    %.4f\n', mu_hat);

% --- 5.3 绘图对比 ---
figure('Name', 'Empirical Analysis', 'Color', 'w');
plot(real_time, real_prices, 'k-', 'LineWidth', 2); hold on;

% 用估计的参数模拟 20 条预测轨道
S_forecast = zeros(real_data_len, 20);
S_forecast(1,:) = real_prices(1);
dW_forecast = sqrt(dt_real) * randn(real_data_len-1, 20);
W_forecast = [zeros(1, 20); cumsum(dW_forecast)];
time_mat_real = repmat(real_time, 1, 20);

% 预测模型
S_forecast = real_prices(1) .* exp((mu_hat - 0.5*sigma_hat^2).*time_mat_real + sigma_hat.*W_forecast);

plot(real_time, S_forecast, 'r', 'Color', [1 0 0 0.7]); % 半透明红色线
xlabel('Time (Year)');
ylabel('Price');
title('真实股价 vs GBM模型模拟 (基于估计参数)');
legend('Historical Data', 'GBM Simulations', 'Location', 'best');
grid on;