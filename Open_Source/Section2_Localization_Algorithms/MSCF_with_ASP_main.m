%% MSCF-IRLS UWB Localization algorithm
% -------------------------------------------------------------------------
% Statement:
%   Demonstration script for UWB TDoA localization with concentrated anchors,
%   including adaptive anchor subset planning and multi-subset consistency fusion
%
% Paper correspondence:
%   1) Distance and ToF model: Eq. (8)-Eq. (9).
%   2) TDoA observation and LOS UWB error model: Eq. (10)-Eq. (18).
%   3) LinHPS linear model and solution: Eq. (19)-Eq. (24).
%   4) Adaptive subset planning: Eq. (41)-Eq. (52).
%   5) Subset-level weighted LinHPS localization: Eq. (53)-Eq. (55).
%   6) Multi-subset consistency fusion: Eq. (56)-Eq. (59).
%   7) IRLS-Huber robust refinement: Eq. (60)-Eq. (67).
%   8) Circular external target test setting: Eq. (74).
%
% Main tunable parameters:
%   Anchor_class_Num  - number of planned anchor subsets M.
%   SS_Anchor_Num     - number of anchors in each subset Ns.
%   Test_R            - target radius around the anchor-cluster COG.
%   Ment_Num          - number of Monte Carlo trials.
%   q1, q2, b1, b2    - LOS ranging-error coefficients, see Eq. (17).
%   n_sync, a_sync    - optional clock-bias scale and fluctuation ratio.
%   max_iter, tol     - maximum iterations and convergence tolerance for IRLS.
%   delta_h lower bound - robust Huber threshold lower bound, see Eq. (63).
%
% Notes for open-source release:
%   This version only enhances comments and documentation. The executable
%   statements of the original script are intentionally kept unchanged.
%   The script expects Anchor.mat to be available in the MATLAB path or the
%   current working directory.
% -------------------------------------------------------------------------
clc;clear;close all;

%% Load anchor layout and build adaptive subsets
% Anchor.mat should contain a 3-by-N matrix named Anchor.
% Each column is one anchor coordinate [x; y; z], consistent with Eq. (3)-Eq. (4).
% In the manuscript, the optimized 12-anchor layout is listed in Table 1.
load("Anchor.mat");
Anchor_all = Anchor;

% Number of anchor subsets M used by MSCF.
% This corresponds to M in Eq. (42) and the parameter setting in Section 4.4.1.
% Increase this value to use more subset estimates; decrease it for lower cost.
Anchor_class_Num = 5;

% Number of anchors in each subset Ns.
% For 3-D LinHPS, Ns should be at least 4; the paper uses Ns = 7.
SS_Anchor_Num = 7;
Anchor_class = AnchorClassPlan_LinHPS(Anchor_all, Anchor_class_Num, SS_Anchor_Num);

%% Generate circular external target positions
% The following angle grid generates 36 azimuth directions around the COG,
% corresponding to the target-location configuration in Eq. (74).
seita = 0 : 1/18*pi : 35/18*pi;
COG = [0.8;0.56;1.72];  % Designed COG(Consistent with the Region_Center proposed in anchor deployment method)

% Target radius R in Eq. (74). In the paper, R is varied from 5 m to 15 m
% for the RMSE/STD curves. Here the script evaluates one selected radius.
Test_R = 15;
for i = 1:numel(seita)
    Test_MS(:,i) = [COG(1)+Test_R*cos(seita(i));COG(2)+Test_R*sin(seita(i));1];
end

for seita_iter = 1:numel(seita)

%% Monte Carlo simulation for the current azimuth direction
% Ment_Num controls the number of repeated noisy trials for each target point.
% The paper uses 500 independent Monte Carlo simulations in Section 4.4.1.
Ment_Num = 500;

for Ment_iter = 1:Ment_Num
tic

%% Physical constants and LOS ranging-error parameters
% c is the propagation speed used by the ToF model in Eq. (9).
% q1, q2, b1, and b2 are the LOS UWB ranging-error coefficients in Eq. (17).
% Adjust these values if a different UWB chip or calibrated ranging model is used.
c = 3e+8;
q1 = 0.0042;
q2 = 0.01;
b1 = -0.0003;
b2 = 0.0302;

Anchor_class_Num = numel(Anchor_class);

plot3(Anchor_all(1,:),Anchor_all(2,:),Anchor_all(3,:),LineStyle="none",Marker=".",MarkerSize=15)
[~,Anchor_all_Num] = size(Anchor_all);
COG = sum(Anchor_all,2)./Anchor_all_Num;  % The real COG(Real Anchor center)
MS = Test_MS(:,seita_iter);  %Mobile station

%Real distances between Anchors and MS
for i = 1:Anchor_all_Num
    R_real(i) = norm(Anchor_all(:,i)-MS,2);
end

%Real TOFs between Anchors and MS
for i = 1:Anchor_all_Num
    TOF_real(i) = R_real(i)/c;
end

%% Optional clock-bias model parameters
% n_sync and a_sync are used for the optional clock-asynchrony experiment.
% In the current executable line, sync_err is generated but not added because
% the clock-bias line below is commented out. This preserves the LOS baseline.
n_sync = 1e-8;
a_sync = 0.15;
for i = 1:Anchor_all_Num
    sync_err = n_sync*(1-a_sync)+rand(1)*2*n_sync*a_sync;
    u(i) = q1*(R_real(i)/c)+q2/c;  %without clock bias
    % u(i) = q1*(R_real(i)/c)+q2/c+sync_err;  %Activate clock bias
    sigma(i) = b1*(R_real(i)/c)+b2/c;
end

%Measured TOFs between Anchors and MS
for i = 1:Anchor_all_Num
    Noise(i) = normrnd(u(i),sigma(i),1);
    TOF_meas(i) = TOF_real(i)+1*Noise(i);
end

%Measured distances between Anchors and MS
for i = 1:Anchor_all_Num
    R_meas(i) = TOF_meas(i)*c;
end

%Measured TDoAs between Anchors and MS
for i = 1:Anchor_all_Num
    delta_TOF_meas(i) = TOF_meas(i)-TOF_meas(1);
end

%% Subset-level LinHPS localization
% Each planned subset is solved independently. The first anchor in each subset
% is treated as the reference anchor after reference-order optimization.
% This block corresponds to Eq. (53)-Eq. (55).
for class_iter = 1:Anchor_class_Num
    SS_index = Anchor_class{1,class_iter};
    SS_Anchor = zeros(3,1);
    SS_TOF_meas = zeros(1);
    SS_R_meas = zeros(1);
    SS_delta_TOF_meas = zeros(1);

    SS_Anchor(:,1) = Anchor_all(:,SS_index(1));
    SS_TOF_meas(1) = TOF_meas(SS_index(1));
    for SS_iter = 2:numel(SS_index)
        SS_Anchor = [SS_Anchor,Anchor_all(:,SS_index(SS_iter))];
        SS_TOF_meas = [SS_TOF_meas,TOF_meas(SS_index(SS_iter))];
    end
    
    for i = 1:numel(SS_index)
        SS_delta_TOF_meas(i) = SS_TOF_meas(i)-SS_TOF_meas(1);
    end
    [~,SS_Anchor_Num] = size(SS_Anchor);
    
    sigma_SS = zeros(1);
    for i = 1:SS_Anchor_Num
        SS_R_meas(i) = SS_TOF_meas(i)*c;
    end

    for i = 1:SS_Anchor_Num
        sigma_SS(i) = b1*(SS_R_meas(i)/c)+b2/c;
    end

% Construct the LinHPS row vector e_ij and scalar b_ij.
% These expressions correspond to Eq. (19)-Eq. (21), and then form E*s = b
% in Eq. (22). The final weighted LS/pseudoinverse solution follows Eq. (23)-Eq. (24)
% and the subset-level weighted form in Eq. (55).
    e = cell(SS_Anchor_Num,SS_Anchor_Num);
    for i = 1:SS_Anchor_Num
        for j = 1:SS_Anchor_Num
            e{i,j} = 2*c*(SS_delta_TOF_meas(j)*(SS_Anchor(:,i)-SS_Anchor(:,1)) ...
                -SS_delta_TOF_meas(i)*(SS_Anchor(:,j)-SS_Anchor(:,1)));
        end
    end

    b = cell(SS_Anchor_Num,SS_Anchor_Num);
    for i = 1:SS_Anchor_Num
        for j = 1:SS_Anchor_Num
            b{i,j} = c*(SS_delta_TOF_meas(i)*(c^2*SS_delta_TOF_meas(j)^2-norm(SS_Anchor(:,j),2)^2) ...
                +(SS_delta_TOF_meas(i)-SS_delta_TOF_meas(j))*norm(SS_Anchor(:,1),2)^2 ...
                +SS_delta_TOF_meas(j)*(norm(SS_Anchor(:,i),2)^2-c^2*SS_delta_TOF_meas(i)^2));
        end
    end

    k = 1;
    E = zeros(1,3);
    for i = 2 : SS_Anchor_Num-1
        for j = i+1 : SS_Anchor_Num
            E(k,:) = e{i,j}';
            sigma_SS_matrix(k,k) = min(max(1/max((c*SS_delta_TOF_meas(j))^2 ...
                *(sigma_SS(i)^2+sigma_SS(1)^2)+ ...
                (c*SS_delta_TOF_meas(i))^2*(sigma_SS(j)^2+ ...
                sigma_SS(1)^2), 1e-12), 1e-3), 1e3);
            k = k+1;
        end
    end

    k = 1;
    B = zeros(1);
    for i = 2 : SS_Anchor_Num-1
        for j = i+1 : SS_Anchor_Num
            B(k,:) = b{i,j};
            k = k+1;
        end
    end

    if rcond(E'*sigma_SS_matrix*E)<1e-14
        S_meas_class(:,class_iter) = pinv(E)*B;
    else
        S_meas_class(:,class_iter) = inv(E'*sigma_SS_matrix*E)*E'*sigma_SS_matrix*B;
    end
end

%% Multi-subset consistency fusion
% The subset estimates are fused according to their geometric consistency with
% the global measured range information. This block implements Eq. (56)-Eq. (59).
% Estimates with smaller consistency deviation receive larger weights.
S_meas_avg = sum(S_meas_class,2)/Anchor_class_Num;
for devia_i = 1:Anchor_class_Num
    devia_dist(devia_i) = abs(norm(S_meas_class(:,devia_i)-COG,2)-sum(R_meas)/Anchor_all_Num);
end

delta_d = devia_dist;
k_sigmoid = 1;
for w_iter = 1:Anchor_class_Num
    delta_d(w_iter) = 1/delta_d(w_iter);
end

for w_iter = 1:Anchor_class_Num
    w(w_iter) = delta_d(w_iter)/sum(delta_d);
end

for w_iter = 1:Anchor_class_Num
    S_meas_class(:,w_iter) = w(w_iter)*S_meas_class(:,w_iter);
end
S_meas = sum(S_meas_class,2);

%% IRLS-Huber robust refinement
% The fused result is used as the initial value s0 for robust fine localization.
% The refinement uses all anchors and measured ranges, corresponding to Eq. (60)-Eq. (67).
% S_meas = refine_IRLS_Huber(S_meas, Anchor_all, R_meas);

error = norm(S_meas(1:3,:)-MS,2);
error_Ment(Ment_iter) = error;
S_meas_Ment(1:3,Ment_iter) = S_meas(1:3,:);
runtime_Ment(Ment_iter) = toc;
end

S_meas_avg = sum(S_meas_Ment,2)/Ment_Num;
error_Ment_avg = sum(error_Ment)/Ment_Num;
sigma_error = sqrt(sum((error_Ment(1,:)-error_Ment_avg).^2,2)/(Ment_Num-1));
runtime_Ment_avg = sum(runtime_Ment)/Ment_Num;

% Store the outputs for the current azimuth angle. These quantities are used
% to compute and visualize the evaluation metrics in Section 4 of the paper.
%
% LinHPS_CRLB(:,seita_iter)      : Cramér-Rao lower bound, corresponding to
%                                 the theoretical reference metric in Eq. (69).
% S_meas_seita(:,seita_iter)     : averaged estimated target position over
%                                 repeated Monte Carlo trials.
% error_avg_seita(seita_iter)    : average localization error over Monte Carlo
%                                 trials, used for the RMSE analysis in Eq. (68).
% sigma_error_seita(seita_iter)  : standard deviation of the localization error,
%                                 corresponding to the STD metric in Eq. (73).
% runtime_avg_seita(seita_iter)  : average runtime of the online localization
%                                 procedure, used for the runtime comparison.
LinHPS_CRLB(:,seita_iter) = CRLB(Anchor_all,MS);
S_meas_seita(:,seita_iter) = S_meas_avg;
error_avg_seita(seita_iter) = error_Ment_avg;
sigma_error_seita(seita_iter) = sigma_error;
runtime_avg_seita(seita_iter) = runtime_Ment_avg;
end

% -------------------------------------------------------------------------
% Local function: IRLS-Huber robust refinement
% -------------------------------------------------------------------------
% This function refines the MSCF coarse estimate using all range measurements.
% Paper mapping:
%   Residual definition: Eq. (60)
%   Huber objective/loss: Eq. (61)-Eq. (62)
%   MAD-based threshold: Eq. (63)
%   Jacobian: Eq. (64)
%   Huber weights: Eq. (65)
%   Weighted least-squares update: Eq. (66)-Eq. (67)
% Tunable parameters inside this function:
%   max_iter - maximum IRLS iterations, set to 10 in Section 4.4.1.
%   tol      - convergence tolerance epsilon, set to 1e-6 in Section 4.4.1.
%   0.01     - lower bound delta_min for the Huber threshold in Eq. (63).
function S_opt = refine_IRLS_Huber(S_init, Anchor_all, R_meas)
    N = numel(R_meas);
    x = S_init(:);
    
    % Maximum number of IRLS iterations Tmax in Eq. (67).
    % Increasing it may improve convergence in difficult cases but increases runtime.
    max_iter = 10;

    % Convergence threshold epsilon in Eq. (67).
    % A smaller value gives a stricter stopping condition.
    tol = 1e-6;
    
    res0 = zeros(N, 1);
    for i = 1:N
        res0(i) = norm(x - Anchor_all(:,i)) - R_meas(i);
    end
    mad_val = median(abs(res0 - median(res0)));

    % Huber threshold delta in Eq. (63).
    % 1.4826 converts MAD to a Gaussian-consistent scale estimate, and 1.345
    % is the common Huber efficiency constant. The 0.01 term is delta_min.
    delta_h = max(1.345 * 1.4826 * mad_val, 0.01); % 防止过小
    
    for iter = 1:max_iter
        r = zeros(N, 1);
        J = zeros(N, 3);
        
        for i = 1:N
            diff_vec = x - Anchor_all(:,i);
            d = norm(diff_vec);
            if d < 1e-12
                d = 1e-12;
            end
            r(i) = d - R_meas(i);
            J(i,:) = (diff_vec / d)';
        end
       
        abs_r = abs(r);
        w = ones(N, 1);
        idx = abs_r > delta_h;
        w(idx) = delta_h ./ abs_r(idx);
        
        W = diag(w);
        
        JtWJ = J' * W * J;
        if rcond(JtWJ) < 1e-14
            dx = pinv(JtWJ) * (J' * W * (-r));
        else
            dx = JtWJ \ (J' * W * (-r));
        end
        
        x = x + dx;
        
        if x(3) < 0
            x(3) = 0;
        end
        
        if norm(dx) < tol
            break;
        end
    end
    
    S_opt = x;
end

% -------------------------------------------------------------------------
% Local function: LinHPS-oriented adaptive anchor subset planning
% -------------------------------------------------------------------------
% This function generates M anchor subsets with Ns anchors per subset.
% Paper mapping:
%   Complete anchor index set and subset definition: Eq. (41)-Eq. (42)
%   Seed selection with reuse penalty: Eq. (43)-Eq. (44)
%   Candidate subset expansion: Eq. (45)-Eq. (48)
%   Geometric cost and planning cost: Eq. (49)-Eq. (52)
% Main tunable weights in this implementation:
%   0.3  - reuse penalty in seed_score, corresponding to lambda_u.
%   5e2  - center-offset penalty weight, corresponding to lambda_1.
%   5e3  - reuse penalty weight, corresponding to lambda_2.
%   5e3  - diversity penalty weight, corresponding to lambda_3.
function Anchor_class = AnchorClassPlan_LinHPS(Anchor_all, Anchor_class_Num, SS_Anchor_Num)

[row_num, Anchor_all_Num] = size(Anchor_all);
if row_num ~= 3
    % error('Anchor_all 应为 3*N');
    error('Anchor_all should be 3*N');
end
if Anchor_class_Num < 1 || Anchor_class_Num ~= floor(Anchor_class_Num)
    % error('Anchor_class_Num 必须为正整数');
    error('Anchor_class_Num should be a positive integer');
end
if SS_Anchor_Num < 4 || SS_Anchor_Num ~= floor(SS_Anchor_Num)
    % error('SS_Anchor_Num 必须为不小于4的正整数');
    error('SS_Anchor_Num should be a positive integer not less than 4');
end
if Anchor_all_Num < SS_Anchor_Num
    % error('总基站数不足');
    error('insufficient total anchors');
end
if any(Anchor_all(3,:) < 0)
    % error('基站高度不能为负');
    error('The height of the anchors cannot be negative');
end

COG = mean(Anchor_all, 2);
MS_test_set = GenerateOuterMSPoints(Anchor_all);

Anchor_class = cell(1, Anchor_class_Num);
use_count = zeros(1, Anchor_all_Num);

for class_iter = 1:Anchor_class_Num
    seed_score = zeros(1, Anchor_all_Num);
    for i = 1:Anchor_all_Num
        seed_score(i) = norm(Anchor_all(:,i)-COG,2) - 0.3*use_count(i);
    end
    [~, seed_idx] = max(seed_score);

    SS_index = seed_idx;

    while numel(SS_index) < SS_Anchor_Num
        candidate_set = setdiff(1:Anchor_all_Num, SS_index);

        best_cost = inf;
        best_idx = candidate_set(1);

        for k = 1:numel(candidate_set)
            idx_try = [SS_index, candidate_set(k)];
            idx_try = BestRefOrder_LinHPS(idx_try, Anchor_all, MS_test_set);

            % 1) cost_1_LinHPS
            cond_cost = SubsetLinHPSCost(idx_try, Anchor_all, MS_test_set);

            % 2) cost_2_bias_of_COG
            COG_SS = mean(Anchor_all(:,idx_try), 2);
            center_cost = norm(COG_SS - COG, 2);

            % 3) cost_3_reuse
            reuse_cost = use_count(candidate_set(k));

            % 4) cost_4_similarity
            diversity_cost = 0;
            if class_iter > 1
                for m = 1:class_iter-1
                    overlap_ratio = numel(intersect(idx_try, Anchor_class{1,m})) / numel(idx_try);
                    diversity_cost = diversity_cost + overlap_ratio;
                end
                diversity_cost = diversity_cost / (class_iter-1);
            end

            total_cost = cond_cost ...
                       + 5e2 * center_cost ...
                       + 5e3 * reuse_cost ...
                       + 5e3 * diversity_cost;

            if total_cost < best_cost
                best_cost = total_cost;
                best_idx = candidate_set(k);
            end
        end

        SS_index = [SS_index, best_idx];
    end

    SS_index = BestRefOrder_LinHPS(SS_index, Anchor_all, MS_test_set);
    Anchor_class{1, class_iter} = SS_index;
    use_count(SS_index) = use_count(SS_index) + 1;
end

end

% -------------------------------------------------------------------------
% Local function: subset geometric stability cost
% -------------------------------------------------------------------------
% This function evaluates the numerical conditioning of a candidate subset by
% building the LinHPS matrix E at representative external test points.
% It implements the geometric-cost idea in Eq. (49): average performance is
% emphasized while the worst case is also considered.
function cond_cost = SubsetLinHPSCost(SS_index, Anchor_all, MS_test_set)

c = 3e8;
test_num = size(MS_test_set, 2);
cond_list = zeros(1, test_num);

Anchor = Anchor_all(:, SS_index);
[~, Anchor_Num] = size(Anchor);

if Anchor_Num < 4
    cond_cost = 1e12;
    return;
end

for k = 1:test_num
    MS = MS_test_set(:,k);

    TOF = zeros(1, Anchor_Num);
    for i = 1:Anchor_Num
        TOF(i) = norm(Anchor(:,i) - MS, 2) / c;
    end
    delta_TOF = TOF - TOF(1);

    row_num = (Anchor_Num-1)*(Anchor_Num-2)/2;
    E = zeros(row_num, 3);

    kk = 1;
    for i = 2:Anchor_Num-1
        for j = i+1:Anchor_Num
            eij = 2*c*(delta_TOF(j)*(Anchor(:,i)-Anchor(:,1)) ...
                - delta_TOF(i)*(Anchor(:,j)-Anchor(:,1)));
            E(kk,:) = eij';
            kk = kk + 1;
        end
    end

    G = E' * E;
    if rank(E) < 3 || rcond(G) < 1e-12
        cond_list(k) = 1e12;
    else
        cond_list(k) = cond(G);   
    end
end

% Geometric cost weighting alpha in Eq. (49).
% Here alpha = 0.9 emphasizes average stability and 0.1 keeps worst-case awareness.
cond_cost = 0.9 * mean(cond_list) + 0.1 * max(cond_list);

end

% -------------------------------------------------------------------------
% Local function: representative external test-point generation
% -------------------------------------------------------------------------
% The generated points are used only for subset planning, not for final testing.
% They approximate external target directions around the concentrated anchor cluster,
% supporting the ASP geometric evaluation described in Section 3.2.2.
function MS_test_set = GenerateOuterMSPoints(Anchor_all)

x_min = min(Anchor_all(1,:));
x_max = max(Anchor_all(1,:));
y_min = min(Anchor_all(2,:));
y_max = max(Anchor_all(2,:));
z_max = max(Anchor_all(3,:));

COG = mean(Anchor_all, 2);

Lx = x_max - x_min;
Ly = y_max - y_min;
% Inner/outer radii for representative ASP test points.
% These values can be adjusted when the expected target region is closer or farther
% from the anchor cluster.
r1 = 0.5 * max(Lx, Ly);
r2 = 10.0 * max(Lx, Ly);

% Azimuth sampling for ASP evaluation. Increasing the number of angles makes
% subset planning more comprehensive but also more computationally expensive.
theta_vec = linspace(0, 2*pi, 12);
theta_vec(end) = [];
z_vec = [0.5, 1.5, min(3.0, max(0.5, z_max+0.3))];

MS_test_set = [];
for r = [r1, r2]
    for iz = 1:numel(z_vec)
        for it = 1:numel(theta_vec)
            x = COG(1) + r*cos(theta_vec(it));
            y = COG(2) + r*sin(theta_vec(it));
            z = z_vec(iz);
            MS_test_set = [MS_test_set, [x;y;z]];
        end
    end
end

end

% -------------------------------------------------------------------------
% Local function: reference-anchor ordering optimization
% -------------------------------------------------------------------------
% LinHPS is sensitive to the reference anchor used to construct TDoA values.
% This function tries each anchor in the subset as the reference anchor and keeps
% the ordering with the smallest subset geometric cost, as described in Section 3.2.3.
function SS_index_best = BestRefOrder_LinHPS(SS_index, Anchor_all, MS_test_set)

best_cost = inf;
SS_index_best = SS_index;

for i = 1:numel(SS_index)
    idx_try = SS_index;
    idx_try([1,i]) = idx_try([i,1]);

    cost_try = SubsetLinHPSCost(idx_try, Anchor_all, MS_test_set);

    if cost_try < best_cost
        best_cost = cost_try;
        SS_index_best = idx_try;
    end
end

end

% -------------------------------------------------------------------------
% Local function: CRLB calculation
% -------------------------------------------------------------------------
function LinHPS_CRLB = CRLB(Anchor,MS)

[~,Anchor_Num] = size(Anchor);
b1 = -0.0003;
b2 = 0.0302;

for i = 1:Anchor_Num
    R_real(i) = norm(Anchor(:,i)-MS,2);
end
sum_R_real = sum(R_real);
sigma_CRLB = sqrt(2)*(b1*(sum_R_real/Anchor_Num)+b2);

for J_row_iter = 1:3
    for J_col_iter = 1:3
        sum_J = 0;
        for i = 1 : Anchor_Num-1
            for j = i:Anchor_Num
                sum_J = sum_J+((MS(J_row_iter)-Anchor(J_row_iter,i))/R_real(i)- ...
                    (MS(J_row_iter)-Anchor(J_row_iter,j))/R_real(j))* ...
                    ((MS(J_col_iter)-Anchor(J_col_iter,i))/R_real(i)- ...
                    (MS(J_col_iter)-Anchor(J_col_iter,j))/R_real(j));
            end
        end
        J_matrix(J_row_iter,J_col_iter) = sum_J/(sigma_CRLB^2);
    end
end
LinHPS_CRLB = sqrt(trace(inv(J_matrix)));

end


