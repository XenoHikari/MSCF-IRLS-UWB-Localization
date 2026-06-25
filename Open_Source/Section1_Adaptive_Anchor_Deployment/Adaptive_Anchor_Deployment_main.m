clc;
clear;
close all;

% =========================================================================
% Adaptive Anchor Deployment for Concentrated UWB-TDoA Localization
% =========================================================================
% This script implements the task-driven dense anchor deployment experiment
% used for the concentrated-anchor external-localization scenario.
%
% Code role in the paper:
%   - This file corresponds to the adaptive anchor deployment module in
%     Section 3.1 of the manuscript.
%   - The optimized anchor layout is evaluated by the LinHPS localization
%     error, centroid offset, and geometric condition number.
%   - The comprehensive deployment metric corresponds to the structure of
%     Eq. (35) - Eq. (40).
%
% Main outputs:
%   Anchor      : optimized 3-D anchor coordinates, size 3 x Anchor_Num.
%   best_metric : minimum comprehensive deployment evaluation value.
%   Record      : optimization history, including mean localization error,
%                 centroid offset penalty, condition-number penalty, and
%                 total metric value.
%
% Open-source note:
%   The implementation is written as a self-contained MATLAB script with
%   local functions. It can be used to reproduce the anchor deployment
%   optimization experiment or to generate anchor layouts for subsequent
%   LinHPS / MSCF-IRLS simulations.
%
% Adjustable parameters:
%   Length, Width, Height : physical size of the concentrated deployment box.
%   Anchor_Num            : number of anchors in the cluster.
%   Region_Center         : desired center of the anchor deployment region.
%   Outer_Range           : half-size of the external target activity region.
%   Search_Num            : number of candidate layouts to be evaluated.
%   MS_Num                : number of external test points per candidate.
%   lambda_cog            : weight of centroid-offset penalty.
%   lambda_cond           : weight of geometric condition-number penalty.
% =========================================================================

% save("Anchor.mat","Anchor");
% rng('shuffle');

% -------------------------------------------------------------------------
% Scenario configuration.
% The anchor cluster is constrained in a compact cuboid deployment region,
% corresponding to the concentrated anchor model in Eq. (1). The center of
% this region is denoted by c in Eq. (2). The external target activity
% region is used to generate random MS test points outside the anchor group.
% -------------------------------------------------------------------------

% Deployment region: 1.5m × 1.5m × 4m
Length = 1.5;
Width = 1.5;
Height = 4.0;

Anchor_Num = 12;

% Designed anchor center
Region_Center = [0.8; 0.56; 1.72];

% MS generated region
Outer_Range = [18; 18; 3];

% -------------------------------------------------------------------------
% Search for the best dense anchor layout.
% The function AutoAnchorDeployDense_Opt generates multiple candidate anchor
% clusters and evaluates each one using:
%   1) average LinHPS localization error, corresponding to Eq. (35);
%   2) centroid offset penalty, corresponding to Eq. (36);
%   3) geometric condition-number penalty, corresponding to Eq. (37)-Eq. (39).
% The final metric follows the weighted form of Eq. (40).
% -------------------------------------------------------------------------

[Anchor, best_metric, Record] = AutoAnchorDeployDense_Opt(Length, Width, Height, Anchor_Num, Region_Center, Outer_Range);
% disp('====== 最优Anchor布设结果 ======');
disp('====== Optimal Anchors ======');
disp(Anchor);
% fprintf('最优综合评价指标 best_metric = %.6f\n', best_metric);
fprintf('best_metric = %.6f\n', best_metric);

% -------------------------------------------------------------------------
% A single external MS point is generated to provide a direct visual and
% numerical check of the optimized anchor layout. The LinHPS function below
% uses the TDoA observation model and the linear hyperbolic equation system
% to estimate the target position.
% -------------------------------------------------------------------------
MS = GenerateMSOutsideAnchorGroup(Region_Center, Length, Width, Height, Outer_Range);
[S_meas, COG] = LinHPS(Anchor, MS);
pos_err = norm(S_meas - MS, 2);

% disp('====== 单次测试结果 ======');
disp('====== Single target test result ======');
% disp('真实MS = ');
disp('True MS = ');
disp(MS);
% disp('估计S_meas = ');
disp('Estimated S_meas = ');
disp(S_meas);
disp('COG = ');
disp(COG);
% fprintf('单次定位误差 = %.6f m\n', pos_err);
fprintf('Single test error = %.6f m\n', pos_err);

% -------------------------------------------------------------------------
% Monte Carlo validation over randomly generated external target positions.
% The error statistics provide a preliminary assessment of whether the
% optimized anchor layout is suitable for external localization. Invalid
% numerical solutions are excluded from the final statistics.
% -------------------------------------------------------------------------
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

% fprintf('\n====== Monte Carlo统计结果 ======\n');
% fprintf('有效测试次数 = %d / %d\n', length(Err_valid), Test_Num);
% fprintf('平均误差 = %.6f m\n', mean(Err_valid));
% fprintf('RMSE = %.6f m\n', sqrt(mean(Err_valid.^2)));
% fprintf('最大误差 = %.6f m\n', max(Err_valid));
% fprintf('最小误差 = %.6f m\n', min(Err_valid));

fprintf('\n====== Monte Carlo result ======\n');
fprintf('Test number = %d / %d\n', length(Err_valid), Test_Num);
fprintf('average error = %.6f m\n', mean(Err_valid));
fprintf('RMSE = %.6f m\n', sqrt(mean(Err_valid.^2)));
fprintf('max error = %.6f m\n', max(Err_valid));
fprintf('min error = %.6f m\n', min(Err_valid));

% -------------------------------------------------------------------------
% Check the actual minimum inter-anchor spacing of the optimized layout.
% This verifies whether the generated dense deployment satisfies the spacing
% constraint used during candidate generation.
% -------------------------------------------------------------------------
d_min_real = inf;
for i = 1:Anchor_Num-1
    for j = i+1:Anchor_Num
        d_ij = norm(Anchor(:,i) - Anchor(:,j), 2);
        if d_ij < d_min_real
            d_min_real = d_ij;
        end
    end
end
% fprintf('实际最小基站间距 = %.6f m\n', d_min_real);
fprintf('minimal anchor distances = %.6f m\n', d_min_real);

% -------------------------------------------------------------------------
% Visualize the optimized anchor cluster, the true MS position, the estimated
% MS position, the centroid of the anchor group, and the deployment cuboid.
% This figure corresponds to the deployment visualization used to explain
% the spatial structure of the optimized anchors.
% -------------------------------------------------------------------------
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
title('Optimized concentrated Anchor Deployment for LinHPS');

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

% -------------------------------------------------------------------------
% Plot the total deployment metric for all candidate anchor layouts.
% A smaller value indicates better combined localization performance,
% centroid consistency, and geometric stability.
% -------------------------------------------------------------------------
figure;
plot(Record.metric, 'LineWidth', 1.5);
grid on;
% xlabel('候选Anchor编号');
xlabel('Candidate Anchor number');
% ylabel('综合评价指标');
ylabel('Evaluation index');
title('Anchor Optimization Process');

% -------------------------------------------------------------------------
% Plot the Monte Carlo localization-error histogram under the optimized
% anchor layout. This is useful for checking whether large-error samples
% occupy only a small proportion of all random tests.
% -------------------------------------------------------------------------
figure;
histogram(Err_valid, 20);
grid on;
% xlabel('定位误差 / m');
xlabel('Localization error / m');
% ylabel('频数');
ylabel('Frequency');
title('Monte Carlo Localization Error Histogram');

% =========================================================================
% Local functions
% =========================================================================
% The following functions are placed in the same script for easier
% open-source reproduction. They include:
%   1) AutoAnchorDeployDense_Opt      : task-driven anchor deployment search.
%   2) GenerateMSOutsideAnchorGroup   : external target-point generation.
%   3) LinHPS                         : baseline TDoA LinHPS solver.
%   4) GenerateDenseAnchorCandidate   : candidate anchor-layout generation.
% =========================================================================
function [Anchor_best, best_metric, Record] = AutoAnchorDeployDense_Opt(Length, Width, Height, Anchor_Num, Region_Center, Outer_Range)
if nargin < 6
    % error('输入参数不足，应为 Length, Width, Height, Anchor_Num, Region_Center, Outer_Range');
    error('Insufficient inputs');
end

% ---------------------------------------------------------------------
% Search_Num controls the number of candidate anchor layouts. Increasing
% it may improve the chance of finding a better layout but increases
% offline computation.
%
% MS_Num controls the number of randomly sampled external target points
% used to estimate the average localization error of each candidate
% layout. A larger value gives a more stable estimate of Eq. (35), at the
% cost of longer offline optimization time.
%
% lambda_cog and lambda_cond are scale-balancing coefficients in the
% comprehensive metric of Eq. (40). They are empirical parameters rather
% than theoretically optimal constants.
% ---------------------------------------------------------------------

Search_Num = 100;
MS_Num = 20000;
lambda_cog = 0.15;
lambda_cond = 0.02;

best_metric = inf;
Anchor_best = [];

Record.metric = zeros(1, Search_Num);
Record.mean_err = zeros(1, Search_Num);
Record.cog_err = zeros(1, Search_Num);
Record.cond_penalty = zeros(1, Search_Num);

for n = 1:Search_Num
    % Candidate generation follows the constrained dense-deployment idea:
    % anchors are sampled inside the deployment cuboid while satisfying
    % the minimum spacing constraint similar to Eq. (28).
    Anchor = GenerateDenseAnchorCandidate(Length, Width, Height, Anchor_Num, Region_Center);

    % External MS samples approximate the target test-point set used for
    % deployment evaluation. The mean error serves as the localization
    % performance term J_loc(A) in Eq. (35).
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

    % Penalize excessive deviation between the generated anchor centroid
    % and the desired deployment-region center, corresponding to Eq. (36).
    COG = sum(Anchor, 2) / Anchor_Num;
    cog_err = norm(COG - Region_Center, 2);

    % Construct a centralized anchor matrix and evaluate the condition
    % number of the associated geometry matrix, corresponding to
    % Eq. (37)-Eq. (39). A large condition number indicates possible
    % geometric degeneration or poor numerical stability.
    A_centered = Anchor - mean(Anchor, 2);
    G = A_centered * A_centered.';
    cond_penalty = cond(G);
    if ~isfinite(cond_penalty)
        cond_penalty = 1e6;
    end

    % Comprehensive deployment score, corresponding to Eq. (40):
    % metric = localization error + centroid penalty + conditioning
    % penalty. The candidate with the minimum metric is retained.
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
    % error('优化失败，未找到有效的Anchor布设。');
    error('Optimization failed');
end
end


function MS = GenerateMSOutsideAnchorGroup(Anchor_Center, Anchor_Length, Anchor_Width, Anchor_Height, Outer_Range)
    if nargin < 5
        % error('输入参数不足。');
        error('Insufficient inputs');
    end

    xc = Anchor_Center(1);
    yc = Anchor_Center(2);
    zc = Anchor_Center(3);

    inner_x_min = xc - Anchor_Length/2;
    inner_x_max = xc + Anchor_Length/2;
    inner_y_min = yc - Anchor_Width/2;
    inner_y_max = yc + Anchor_Width/2;
    inner_z_min = max(0, zc - Anchor_Height/2);
    inner_z_max = zc + Anchor_Height/2;

    % The buffer excludes points too close to the anchor cuboid. This ensures
    % that the generated test points represent the external-localization
    % scenario considered in the paper.
    buffer = 0.5;

    outer_x_min = xc - Outer_Range(1);
    outer_x_max = xc + Outer_Range(1);
    outer_y_min = yc - Outer_Range(2);
    outer_y_max = yc + Outer_Range(2);
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

    % error('无法生成位于基站群外侧的MS，请检查Outer_Range设置。');
    error('Generation failed');
end


function [S_meas,COG] = LinHPS(Anchor,MS)
% LinHPS
% Baseline Linear Hyperbolic Positioning System solver.
%
% This function is used only as the basic localizer inside the deployment
% evaluation. For each candidate anchor layout, LinHPS estimates the MS
% position, and the resulting error is accumulated as the task-driven
% deployment objective.
%
% Formula correspondence:
%   R_real       : anchor-to-target distance, corresponding to Eq. (8).
%   TOF_real     : ideal propagation time, corresponding to Eq. (9).
%   delta_TOF    : TDoA observation relative to Anchor 1, Eq. (10)-Eq. (11).
%   q1,q2,b1,b2  : distance-dependent UWB LOS error model, Eq. (15)-Eq. (17).
%   E and B      : LinHPS linear system E*s = b, Eq. (19)-Eq. (22).
%   S_meas       : pseudoinverse solution, corresponding to Eq. (24).
%
% Adjustable parameters:
%   q1, q2, b1, b2 can be changed if a different UWB ranging module or
%   measurement calibration dataset is used.
%
% Speed of light used to convert ToF and TDoA into distance quantities.
c = 3e+8;
% UWB LOS ranging-error coefficients. These values are consistent with the
% distance-dependent measurement model used in the manuscript.
q1 = 0.0042;
q2 = 0.01;
b1 = -0.0003;
b2 = 0.0302;
[~,Anchor_Num] = size(Anchor);
COG = sum(Anchor,2)./Anchor_Num;  % The real COG(Anchor Center)

% True range between each anchor and the MS, corresponding to Eq. (8).
for i = 1:Anchor_Num
    R_real(i) = norm(Anchor(:,i)-MS,2);
end

% Ideal time of flight, corresponding to Eq. (9).
for i = 1:Anchor_Num
    TOF_real(i) = R_real(i)/c;
end

% Distance-dependent LOS measurement bias and standard deviation,
% corresponding to Eq. (15)-Eq. (17).
n_sync = 1e-8;
a_sync = 0.1;
for i = 1:Anchor_Num
    sync_err = n_sync*(1-a_sync)+rand(1)*2*n_sync*a_sync;
    u(i) = q1*(R_real(i)/c)+q2/c;
    %u(i) = q1*(R_real(i)/c)+q2/c+sync_err;
    sigma(i) = b1*(R_real(i)/c)+b2/c;
end

% Measured ToF is obtained by adding the simulated UWB LOS noise to the
% ideal ToF. The measured range R_meas corresponds to Eq. (18).
for i = 1:Anchor_Num
    Noise(i) = normrnd(u(i),sigma(i),1);
    TOF_meas(i) = TOF_real(i)+1*Noise(i);
end
for i = 1:Anchor_Num
    R_meas(i) = TOF_meas(i)*c;
end

% TDoA measurement relative to Anchor 1, corresponding to Eq. (10)-Eq. (11).
for i = 1:Anchor_Num
    delta_TOF_meas(i) = TOF_meas(i)-TOF_meas(1);
end

% Construct the LinHPS linear equation system. The intermediate vector e_ij
% and scalar b_ij correspond to Eq. (19) and Eq. (20), and the stacked system
% corresponds to Eq. (21)-Eq. (22).
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

% Solve the LinHPS system using the Moore-Penrose pseudoinverse. This is
% numerically more tolerant when E is ill-conditioned, corresponding to
% Eq. (24).
S_meas = pinv(E)*B;
end


function Anchor = GenerateDenseAnchorCandidate(Length, Width, Height, Anchor_Num, Region_Center)
if nargin < 5
    % error('输入参数不足，应为 Length, Width, Height, Anchor_Num, Region_Center');
    error('Insufficient inputs');
end

% Estimate the average volumetric spacing from the deployment volume.
% This provides a scale-aware basis for the minimum inter-anchor spacing.
V = Length * Width * Height;
d_avg = (V / Anchor_Num)^(1/3);

% Minimum inter-anchor spacing constraint, corresponding to the role of
% d_min in Eq. (28). The coefficient 0.6 can be adjusted to control the
% compactness of the generated anchor cluster.
d_min = 0.6 * d_avg;

% Boundary margins prevent anchors from being placed too close to the
% deployment-box boundary, corresponding to the practical margin idea in
% Eq. (26)-Eq. (27).
margin_x = min(0.05 * Length, d_min / 2);
margin_y = min(0.05 * Width,  d_min / 2);
margin_z = min(0.05 * Height, d_min / 2);

% Actual feasible deployment region after applying boundary margins.
% The lower bound of z is clipped to keep anchor heights non-negative.
x_min = Region_Center(1) - Length/2 + margin_x;
x_max = Region_Center(1) + Length/2 - margin_x;
y_min = Region_Center(2) - Width/2  + margin_y;
y_max = Region_Center(2) + Width/2  - margin_y;
z_min = max(0 + margin_z, Region_Center(3) - Height/2 + margin_z);
z_max = Region_Center(3) + Height/2 - margin_z;
if z_max <= z_min
    % error('基站布设高度范围无效，请检查Region_Center(3)与Height设置。');
    error('Invalid z_max');
end

% Build a 3-D candidate grid according to the aspect ratio of the
% deployment box. This corresponds to the grid-based candidate generation
% described by Eq. (29)-Eq. (30).
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
% Generate all grid candidate points and randomly shuffle their order so
% that different calls can produce different feasible anchor layouts.
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

% Select anchors from the shuffled grid candidates while enforcing the
% minimum spacing constraint.
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

% If the grid candidates are not sufficient, uniformly sample additional
% random points in the feasible deployment region, similar to Eq. (31).
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
    % error('当前参数下无法生成满足最小间距约束的Anchor，请调整区域尺寸、Anchor_Num或最小间距比例。');
    error('Invalid Anchor_Num or count');
end

% Apply small random perturbations to break overly regular grid symmetry
% while keeping the anchors inside the feasible region and satisfying the
% minimum spacing constraint, corresponding to Eq. (32)-Eq. (33).
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