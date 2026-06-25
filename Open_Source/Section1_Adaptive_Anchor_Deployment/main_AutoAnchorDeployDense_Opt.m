clc;
clear;
close all;

% save("Anchor.mat","Anchor");
% rng('shuffle');

%% 场景参数
% 基站密集分布空间：1.5m × 1.5m × 4m
Length = 1.5;
Width = 1.5;
Height = 4.0;

Anchor_Num = 12;

% 基站群中心，可产生负值坐标
Region_Center = [0.8; 0.56; 1.72];

% MS活动区域半尺寸（MS位于基站群外侧）
% 例如可覆盖到类似[-15; -10; 5]这样的坐标
Outer_Range = [18; 18; 3];
% for i = 1:1000
%     Test_MS(:, i) = GenerateMSOutsideAnchorGroup(Region_Center, Length, Width, Height, Outer_Range);
% end

%% 优化基站布设
[Anchor, best_metric, Record] = AutoAnchorDeployDense_Opt(Length, Width, Height, Anchor_Num, Region_Center, Outer_Range);
disp('====== 最优Anchor布设结果 ======');
disp(Anchor);
fprintf('最优综合评价指标 best_metric = %.6f\n', best_metric);

%% 单次测试
MS = GenerateMSOutsideAnchorGroup(Region_Center, Length, Width, Height, Outer_Range);
[S_meas, COG] = LinHPS(Anchor, MS);
pos_err = norm(S_meas - MS, 2);

disp('====== 单次测试结果 ======');
disp('真实MS = ');
disp(MS);
disp('估计S_meas = ');
disp(S_meas);
disp('COG = ');
disp(COG);
fprintf('单次定位误差 = %.6f m\n', pos_err);

%% Monte Carlo统计测试
Test_Num = 3000;
Err = zeros(1, Test_Num);

for k = 1:Test_Num
    MS_test = GenerateMSOutsideAnchorGroup(Region_Center, Length, Width, Height, Outer_Range);
    try
        [S_test, ~] = LinHPS(Anchor, MS_test);
        Err(k) = norm(S_test - MS_test, 2);
    catch
        Err(k) = NaN;
    end
end

Err_valid = Err(~isnan(Err) & ~isinf(Err));

fprintf('\n====== Monte Carlo统计结果 ======\n');
fprintf('有效测试次数 = %d / %d\n', length(Err_valid), Test_Num);
fprintf('平均误差 = %.6f m\n', mean(Err_valid));
fprintf('RMSE = %.6f m\n', sqrt(mean(Err_valid.^2)));
fprintf('最大误差 = %.6f m\n', max(Err_valid));
fprintf('最小误差 = %.6f m\n', min(Err_valid));

%% 实际最小基站间距
d_min_real = inf;
for i = 1:Anchor_Num-1
    for j = i+1:Anchor_Num
        d_ij = norm(Anchor(:,i) - Anchor(:,j), 2);
        if d_ij < d_min_real
            d_min_real = d_ij;
        end
    end
end
fprintf('实际最小基站间距 = %.6f m\n', d_min_real);

%% 三维显示
figure;
hold on;
grid on;
box on;

plot3(Anchor(1,:), Anchor(2,:), Anchor(3,:), 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot3(MS(1), MS(2), MS(3), 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 14);
plot3(S_meas(1), S_meas(2), S_meas(3), 'gs', ...
    'MarkerFaceColor', 'g', 'MarkerSize', 10);
plot3(COG(1), COG(2), COG(3), 'kd', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 10);

legend('Anchor', 'MS True', 'MS Estimated', 'COG');
xlabel('x / m');
ylabel('y / m');
zlabel('z / m');
title('Optimized Dense Anchor Deployment for LinHPS');

% 绘制基站群空间包围盒
x1 = Region_Center(1) - Length/2;
x2 = Region_Center(1) + Length/2;
y1 = Region_Center(2) - Width/2;
y2 = Region_Center(2) + Width/2;
z1 = Region_Center(3) - Height/2;
z2 = Region_Center(3) + Height/2;

box_pts = [x1 y1 z1;
           x2 y1 z1;
           x2 y2 z1;
           x1 y2 z1;
           x1 y1 z2;
           x2 y1 z2;
           x2 y2 z2;
           x1 y2 z2]';

edges = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;3 7;4 8];
for e = 1:size(edges,1)
    p1 = box_pts(:, edges(e,1));
    p2 = box_pts(:, edges(e,2));
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k--');
end

view(3);

%% 优化过程
figure;
plot(Record.metric, 'LineWidth', 1.5);
grid on;
xlabel('候选Anchor编号');
ylabel('综合评价指标');
title('Anchor Optimization Process');

%% 误差直方图
figure;
histogram(Err_valid, 20);
grid on;
xlabel('定位误差 / m');
ylabel('频数');
title('Monte Carlo Localization Error Histogram');



%% 局部函数
function [Anchor_best, best_metric, Record] = AutoAnchorDeployDense_Opt(Length, Width, Height, Anchor_Num, Region_Center, Outer_Range)
% AutoAnchorDeployDense_Opt
% 基于LinHPS定位误差评价的密集分布UWB基站自动优化布设
%
% 输入：
%   Length, Width, Height - 基站群布设空间尺寸
%   Anchor_Num            - 基站数量
%   Region_Center         - 基站群中心 [xc; yc; zc]
%   Outer_Range           - MS外侧活动区域半尺寸 [Rx; Ry; Rz]
%
% 输出：
%   Anchor_best           - 最优Anchor，3×Anchor_Num
%   best_metric           - 最优综合评价指标
%   Record                - 优化过程记录

    if nargin < 6
        error('输入参数不足，应为 Length, Width, Height, Anchor_Num, Region_Center, Outer_Range');
    end

    % 搜索参数
    Search_Num = 100;     % 候选Anchor组数
    MS_Num = 20000;         % 每组Anchor的随机MS测试点数
    lambda_cog = 0.15;   % COG偏移惩罚
    lambda_cond = 0.02;  % 几何病态惩罚

    best_metric = inf;
    Anchor_best = [];

    Record.metric = zeros(1, Search_Num);
    Record.mean_err = zeros(1, Search_Num);
    Record.cog_err = zeros(1, Search_Num);
    Record.cond_penalty = zeros(1, Search_Num);

    for n = 1:Search_Num
        % 1. 生成一组候选Anchor
        Anchor = GenerateDenseAnchorCandidate(Length, Width, Height, Anchor_Num, Region_Center);

        % 2. 随机生成MS（位于基站群外侧），计算平均误差
        err_sum = 0;
        valid_num = 0;

        for k = 1:MS_Num
            MS = GenerateMSOutsideAnchorGroup(Region_Center, Length, Width, Height, Outer_Range);

            try
                [S_meas, ~] = LinHPS(Anchor, MS);

                if any(isnan(S_meas)) || any(isinf(S_meas))
                    continue;
                end

                pos_err = norm(S_meas - MS, 2);
                if isfinite(pos_err)
                    err_sum = err_sum + pos_err;
                    valid_num = valid_num + 1;
                end
            catch
                continue;
            end
        end

        if valid_num < max(10, round(0.3 * MS_Num))
            mean_err = 1e6;
        else
            mean_err = err_sum / valid_num;
        end

        % 3. COG惩罚项
        COG = sum(Anchor, 2) / Anchor_Num;
        cog_err = norm(COG - Region_Center, 2);

        % 4. 几何条件数惩罚
        A_centered = Anchor - mean(Anchor, 2);
        G = A_centered * A_centered.';
        cond_penalty = cond(G);
        if ~isfinite(cond_penalty)
            cond_penalty = 1e6;
        end

        % 5. 综合评价
        metric = mean_err + lambda_cog * cog_err + lambda_cond * cond_penalty;

        Record.metric(n) = metric;
        Record.mean_err(n) = mean_err;
        Record.cog_err(n) = cog_err;
        Record.cond_penalty(n) = cond_penalty;

        if metric < best_metric
            best_metric = metric;
            Anchor_best = Anchor;
        end
    end

    if isempty(Anchor_best)
        error('优化失败，未找到有效的Anchor布设。');
    end
end


function MS = GenerateMSOutsideAnchorGroup(Anchor_Center, Anchor_Length, Anchor_Width, Anchor_Height, Outer_Range)
% GenerateMSOutsideAnchorGroup
% 在基站群外侧随机生成MS位置，允许负值坐标
%
% 输入：
%   Anchor_Center                        - 基站群中心 [xc; yc; zc]
%   Anchor_Length, Anchor_Width, Anchor_Height - 基站群尺寸
%   Outer_Range                          - MS活动范围半尺寸 [Rx; Ry; Rz]
%
% 输出：
%   MS                                   - 3×1，位于基站群外侧的随机坐标

    if nargin < 5
        error('输入参数不足。');
    end

    xc = Anchor_Center(1);
    yc = Anchor_Center(2);
    zc = Anchor_Center(3);

    inner_x_min = xc - Anchor_Length/2;
    inner_x_max = xc + Anchor_Length/2;
    inner_y_min = yc - Anchor_Width/2;
    inner_y_max = yc + Anchor_Width/2;
    % inner_z_min = zc - Anchor_Height/2;
    inner_z_min = max(0, zc - Anchor_Height/2);
    inner_z_max = zc + Anchor_Height/2;

    % 缓冲带，确保MS在基站群外侧
    buffer = 0.5;

    outer_x_min = xc - Outer_Range(1);
    outer_x_max = xc + Outer_Range(1);
    outer_y_min = yc - Outer_Range(2);
    outer_y_max = yc + Outer_Range(2);
    % outer_z_min = zc - Outer_Range(3);
    outer_z_min = 0;
    outer_z_max = zc + Outer_Range(3);

    max_try = 10000;
    try_num = 0;

    while try_num < max_try
        try_num = try_num + 1;

        MS = [outer_x_min + (outer_x_max - outer_x_min) * rand;
              outer_y_min + (outer_y_max - outer_y_min) * rand;
              outer_z_min + (outer_z_max - outer_z_min) * rand];

        in_inner_box = ...
            (MS(1) >= inner_x_min - buffer) && (MS(1) <= inner_x_max + buffer) && ...
            (MS(2) >= inner_y_min - buffer) && (MS(2) <= inner_y_max + buffer) && ...
            (MS(3) >= inner_z_min - buffer) && (MS(3) <= inner_z_max + buffer);

        if ~in_inner_box
            return;
        end
    end

    error('无法生成位于基站群外侧的MS，请检查Outer_Range设置。');
end


function [S_meas,COG] = LinHPS(Anchor,MS)
c = 3e+8;
q1 = 0.0042;
q2 = 0.01;
b1 = -0.0003;
b2 = 0.0302;
[~,Anchor_Num] = size(Anchor);
COG = sum(Anchor,2)./Anchor_Num;  %COG坐标(Hotspot)

%计算各个Anchor与MS之间的真实距离
for i = 1:Anchor_Num
    R_real(i) = norm(Anchor(:,i)-MS,2);
end
%计算各个Anchor与MS之间的真实信号飞行时间
for i = 1:Anchor_Num
    TOF_real(i) = R_real(i)/c;
end
%LOS情况下噪声情况
%n_sync = 1e-8;
%a_sync = 0.1;
for i = 1:Anchor_Num
    %sync_err = n_sync*(1-a_sync)+rand(1)*2*n_sync*a_sync;
    u(i) = q1*(R_real(i)/c)+q2/c;
    %u(i) = q1*(R_real(i)/c)+q2/c+sync_err;  %时钟不同步的情况下
    sigma(i) = b1*(R_real(i)/c)+b2/c;
end
%计算各个Anchor与MS之间的测量信号飞行时间
for i = 1:Anchor_Num
    %Noise(i) = sigma(i)*randn(1)+u(i);  %噪声Noise符合~N(u,sigma)的正态分布噪声
    Noise(i) = normrnd(u(i),sigma(i),1);
    TOF_meas(i) = TOF_real(i)+1*Noise(i);
end
for i = 1:Anchor_Num
    R_meas(i) = TOF_meas(i)*c;
end
%计算各个Anchor与Anchor_1参考锚节点之间的测量信号飞行时间
for i = 1:Anchor_Num
    delta_TOF_meas(i) = TOF_meas(i)-TOF_meas(1);
end

%Es = b公式 利用元胞矩阵存储矩阵数据
e = cell(Anchor_Num,Anchor_Num);
for i = 1:Anchor_Num
    for j = 1:Anchor_Num
        e{i,j} = 2*c*(delta_TOF_meas(j)*(Anchor(:,i)-Anchor(:,1)) ...
            -delta_TOF_meas(i)*(Anchor(:,j)-Anchor(:,1)));
    end
end

b = cell(Anchor_Num,Anchor_Num);
for i = 1:Anchor_Num
    for j = 1:Anchor_Num
        b{i,j} = c*(delta_TOF_meas(i)*(c^2*delta_TOF_meas(j)^2-norm(Anchor(:,j),2)^2) ...
           +(delta_TOF_meas(i)-delta_TOF_meas(j))*norm(Anchor(:,1),2)^2 ...
           +delta_TOF_meas(j)*(norm(Anchor(:,i),2)^2-c^2*delta_TOF_meas(i)^2));
    end
end

k = 1;
for i = 2 : Anchor_Num-1
    for j = i+1 : Anchor_Num
        E(k,:) = e{i,j}';
        k = k+1;
    end
end

k = 1;
for i = 2 : Anchor_Num-1
    for j = i+1 : Anchor_Num
        B(k,:) = b{i,j};
        k = k+1;
    end
end

S_meas = pinv(E)*B;
end


function Anchor = GenerateDenseAnchorCandidate(Length, Width, Height, Anchor_Num, Region_Center)
% GenerateDenseAnchorCandidate
% 在以Region_Center为中心的长方体区域内生成满足最小间距约束的基站候选解
%
% 输入：
%   Length, Width, Height - 基站群布设空间尺寸
%   Anchor_Num            - 基站数量
%   Region_Center         - 布设空间中心坐标 [xc; yc; zc]
%
% 输出：
%   Anchor                - 3×Anchor_Num，每列为一个基站坐标

    if nargin < 5
        error('输入参数不足，应为 Length, Width, Height, Anchor_Num, Region_Center');
    end

    V = Length * Width * Height;
    d_avg = (V / Anchor_Num)^(1/3);

    % 最小间距约束
    d_min = 0.6 * d_avg;

    % 边界留白
    margin_x = min(0.05 * Length, d_min / 2);
    margin_y = min(0.05 * Width,  d_min / 2);
    margin_z = min(0.05 * Height, d_min / 2);

    % 实际可布设边界
    x_min = Region_Center(1) - Length/2 + margin_x;
    x_max = Region_Center(1) + Length/2 - margin_x;
    y_min = Region_Center(2) - Width/2  + margin_y;
    y_max = Region_Center(2) + Width/2  - margin_y;
    z_min = max(0 + margin_z, Region_Center(3) - Height/2 + margin_z);
    z_max = Region_Center(3) + Height/2 - margin_z;
    if z_max <= z_min
        error('基站布设高度范围无效，请检查Region_Center(3)与Height设置。');
    end

    % 根据空间尺寸比例构造网格
    scale = (Anchor_Num / (Length * Width * Height))^(1/3);
    Nx = max(1, round(Length * scale));
    Ny = max(1, round(Width  * scale));
    Nz = max(1, round(Height * scale));

    while Nx * Ny * Nz < Anchor_Num
        step_x = Length / Nx;
        step_y = Width  / Ny;
        step_z = Height / Nz;

        [~, idx] = max([step_x, step_y, step_z]);
        if idx == 1
            Nx = Nx + 1;
        elseif idx == 2
            Ny = Ny + 1;
        else
            Nz = Nz + 1;
        end
    end

    if Nx == 1
        x_set = Region_Center(1);
    else
        x_set = linspace(x_min, x_max, Nx);
    end

    if Ny == 1
        y_set = Region_Center(2);
    else
        y_set = linspace(y_min, y_max, Ny);
    end

    if Nz == 1
        z_set = Region_Center(3);
    else
        z_set = linspace(z_min, z_max, Nz);
    end

    % 生成网格候选点
    Candidate = [];
    for ix = 1:Nx
        for iy = 1:Ny
            for iz = 1:Nz
                Candidate = [Candidate, [x_set(ix); y_set(iy); z_set(iz)]];
            end
        end
    end

    Candidate = Candidate(:, randperm(size(Candidate, 2)));

    Anchor = zeros(3, Anchor_Num);
    count = 0;

    % 先从规则候选中选
    for i = 1:size(Candidate, 2)
        p = Candidate(:, i);

        if count == 0
            count = count + 1;
            Anchor(:, count) = p;
        else
            dist = vecnorm(Anchor(:, 1:count) - p, 2, 1);
            if all(dist >= d_min)
                count = count + 1;
                Anchor(:, count) = p;
            end
        end

        if count >= Anchor_Num
            break;
        end
    end

    % 不足时随机补点
    max_try = 5000;
    try_num = 0;

    while count < Anchor_Num && try_num < max_try
        try_num = try_num + 1;

        p = [x_min + (x_max - x_min) * rand;
             y_min + (y_max - y_min) * rand;
             z_min + (z_max - z_min) * rand];

        dist = vecnorm(Anchor(:, 1:count) - p, 2, 1);
        if all(dist >= d_min)
            count = count + 1;
            Anchor(:, count) = p;
        end
    end

    if count < Anchor_Num
        error('当前参数下无法生成满足最小间距约束的Anchor，请调整区域尺寸、Anchor_Num或最小间距比例。');
    end

    % 小幅随机抖动
    jitter_x = min((Length / max(Nx,1)) * 0.15, d_min * 0.2);
    jitter_y = min((Width  / max(Ny,1)) * 0.15, d_min * 0.2);
    jitter_z = min((Height / max(Nz,1)) * 0.15, d_min * 0.2);

    for i = 1:Anchor_Num
        p_old = Anchor(:, i);

        for t = 1:50
            p_new = p_old + [2*rand-1; 2*rand-1; 2*rand-1] .* [jitter_x; jitter_y; jitter_z];

            p_new(1) = min(max(p_new(1), x_min), x_max);
            p_new(2) = min(max(p_new(2), y_min), y_max);
            p_new(3) = min(max(p_new(3), z_min), z_max);

            idx_other = [1:i-1, i+1:Anchor_Num];
            if isempty(idx_other)
                Anchor(:, i) = p_new;
                break;
            else
                dist = vecnorm(Anchor(:, idx_other) - p_new, 2, 1);
                if all(dist >= d_min)
                    Anchor(:, i) = p_new;
                    break;
                end
            end
        end
    end
end