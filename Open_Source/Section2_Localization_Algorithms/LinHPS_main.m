% ========================================================================
% LinHPS UWB-TDoA Localization algorithm
% ------------------------------------------------------------------------
% Statement:
%   This script evaluates the baseline Linear Hyperbolic Positioning System
%   (LinHPS) algorithm for 3-D UWB TDoA localization under a concentrated
%   anchor configuration.
%
% Relation to the paper:
%   1) Anchor-target distance model:       Eq. (8)
%   2) ToF model:                          Eq. (9)
%   3) TDoA observation model:             Eq. (10)--Eq. (11)
%   4) LOS distance-dependent noise model: Eq. (13)--Eq. (18)
%   5) LinHPS linear system construction:  Eq. (19)--Eq. (24)
%   6) Circular external target sampling:  Eq. (74)
%   7) RMSE / CRLB / STD metrics:          Eq. (68), Eq. (69), Eq. (73)
%
% Adjustable parameters:
%   Test_R       : horizontal radius of the external target circle.
%   Ment_Num     : number of Monte Carlo trials for each azimuth angle.
%   q1, q2       : mean-bias coefficients of the UWB LOS ToF error model.
%   b1, b2       : standard-deviation coefficients of the UWB LOS ToF error
%                  model.
%   n_sync       : clock-bias fluctuation scale used when the synchronization
%                  error model is enabled.
%   a_sync       : relative fluctuation ratio of the clock-bias term.
%   Anchor       : anchor coordinates loaded from Anchor.mat. Different
%                  subsets can be selected by uncommenting the candidate
%                  anchor-selection lines below.
%
% Notes for open-source release:
%   - This script does not include MSCF, ASP, or IRLS-Huber refinement. It is
%     only the simulation of LinHPS 
% ========================================================================
clc;clear;close all;

%% ====================== Anchor loading and target generation ======================
% Load the anchor coordinates used in the simulation.
% The loaded file should provide the variable Anchor, or related anchor
% variables depending on the user's data file.
%
% Anchor format:
%   Anchor is a 3-by-N matrix, where each column is one anchor coordinate:
%       Anchor(:,i) = [x_i; y_i; z_i]
% This is consistent with the anchor-coordinate definition in Eq. (3)--Eq. (4)
% of the paper.

% the Anchor_all is adopted in reference[19]
% Anchor_all = [  0        0    1.2    1.2   0     0   1.2    1.2    1.2   1.2    1.2    1.2;
%                 0     1.12      0   1.12   0  1.12     0   1.12      0  1.12   0.56   0.56;
%                 0.4    0.4    0.4    0.4   2     2   2.2    2.2    3.6   3.6    3.1    0.4];
% Anchor = Anchor_all;
load("Anchor.mat");

% Generate 36 azimuth angles on a circle around the anchor hotspot center.
% This corresponds to the external circular target generation strategy in
% Eq. (74), where the target height is fixed as 1 m.
seita = 0 : 1/18*pi : 35/18*pi;
COG = [0.8;0.56];  % Designed COG(Consistent with the Region_Center proposed in anchor deployment method)
Test_R = 15;
for i = 1:numel(seita)
    Test_MS(:,i) = [COG(1)+Test_R*cos(seita(i));COG(2)+Test_R*sin(seita(i));1];
end

%% ====================== Azimuth-angle loop ======================
% Each azimuth angle corresponds to one target position on the external
% circle. For each target, repeated Monte Carlo trials are performed to
% evaluate the mean error, STD, CRLB, and average runtime.
for seita_iter = 1:numel(seita)

%% ====================== Monte Carlo simulation loop ======================
% Ment_Num controls the number of repeated noisy measurements for each
% target position. A larger value gives smoother statistical results but
% increases the total simulation time.
Ment_Num = 500;
for Ment_iter = 1:Ment_Num
% Start timing the online localization process for the current Monte Carlo
% trial. The elapsed time is later stored in runtime_Ment.
tic

%% ====================== UWB LOS measurement model parameters ======================
% c is the signal propagation speed, used to convert between ToF and range.
% q1, q2, b1, and b2 correspond to the LOS distance-dependent error model in
% Eq. (15)--Eq. (17). The same parameters are used in the paper simulations.
c = 3e+8;
q1 = 0.0042;
q2 = 0.01;
b1 = -0.0003;
b2 = 0.0302;

%% ====================== Anchor selection ======================
% The first column of Anchor is used as the reference anchor in the LinHPS
% TDoA construction. This corresponds to the reference-anchor setting in
% Eq. (10)--Eq. (11).
%
% The following commented examples can be used to test LinHPS with different
% numbers of selected anchors. For 3-D LinHPS, at least four anchors are
% generally required, while more anchors provide more pairwise equations.


% Visualize the anchor layout. This plot is mainly used for checking the
% concentrated anchor geometry and has no influence on localization results.
plot3(Anchor(1,:),Anchor(2,:),Anchor(3,:),LineStyle="none",Marker=".",MarkerSize=15)
[r_Anchor,Anchor_Num] = size(Anchor);
COG = sum(Anchor,2)./Anchor_Num;  % The real COG(Real Anchor center)
MS = Test_MS(:,seita_iter);

%% ====================== True range and true ToF generation ======================
% The true Euclidean distance between the target and each anchor is computed
% according to Eq. (8). The corresponding true ToF is then computed according
% to Eq. (9).
for i = 1:Anchor_Num
    R_real(i) = norm(Anchor(:,i)-MS,2);
end

for i = 1:Anchor_Num
    TOF_real(i) = R_real(i)/c;
end

%% ====================== LOS ToF noise generation ======================
% The following block implements the LOS UWB measurement error model in
% Eq. (13)--Eq. (18). For each anchor, the ToF error is sampled from a normal
% distribution with distance-dependent mean u(i) and standard deviation
% sigma(i).
%
% Optional clock-bias experiment:
%   The variable sync_err is prepared for non-ideal synchronization tests.
%   To enable the clock-bias setting corresponding to Eq. (75), uncomment the
%   line containing '+sync_err' and comment the current u(i) line.

n_sync = 1e-8;
a_sync = 0.15;
for i = 1:Anchor_Num
    sync_err = n_sync*(1-a_sync)+rand(1)*2*n_sync*a_sync;
    u(i) = q1*(R_real(i)/c)+q2/c;
    %u(i) = q1*(R_real(i)/c)+q2/c+sync_err;  %Activate clock-bias
    sigma(i) = b1*(R_real(i)/c)+b2/c;
end

% Add the sampled ToF noise to the true ToF. The measured range R_meas is
% also computed for metric analysis or for possible extensions.
for i = 1:Anchor_Num
    Noise(i) = normrnd(u(i),sigma(i),1);
    TOF_meas(i) = TOF_real(i)+1*Noise(i);
end

for i = 1:Anchor_Num
    R_meas(i) = TOF_meas(i)*c;
end

%% ====================== TDoA measurement construction ======================
% The measured TDoA is constructed by subtracting the ToF of the first anchor
% from all anchor ToF measurements, as defined in Eq. (10)--Eq. (11).
for i = 1:Anchor_Num
    delta_TOF_meas(i) = TOF_meas(i)-TOF_meas(1);
end

%% ====================== LinHPS linear equation construction ======================
% The following two cell arrays store the intermediate e_ij vector and b_ij
% scalar used by LinHPS. They correspond to Eq. (19) and Eq. (20). The final
% stacked matrix equation is E*s = b, as shown in Eq. (21)--Eq. (22).
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

% Stack e_ij into the coefficient matrix E. Only pairs with i,j >= 2 and
% i < j are used, because anchor 1 has already been selected as the internal
% reference anchor.
k = 1;
for i = 2 : Anchor_Num-1
    for j = i+1 : Anchor_Num
        E(k,:) = e{i,j}';
        k = k+1;
    end
end

% Stack b_ij into the right-hand-side vector B.
k = 1;
for i = 2 : Anchor_Num-1
    for j = i+1 : Anchor_Num
        B(k,:) = b{i,j};
        k = k+1;
    end
end

%% ====================== LinHPS position solution ======================
% Solve the LinHPS linear system. The pseudoinverse form corresponds to
% Eq. (24), and is more numerically stable when E is ill-conditioned or rank
% deficient. The commented normal-equation form corresponds to Eq. (23).
S_meas = pinv(E)*B;
%S_meas = inv(E'*E)*E'*B


%% ====================== Monte Carlo outputs for the current trial ======================
% Compute the Euclidean localization error for the current noisy trial and
% store the estimated position and runtime.
error = norm(S_meas(1:3,:)-MS,2);
error_Ment(Ment_iter) = error;
S_meas_Ment(1:3,Ment_iter) = S_meas(1:3,1);
runtime_Ment(Ment_iter) = toc;
end

%% ====================== Statistics over Monte Carlo trials ======================
% Average the estimated positions and compute the localization-error mean,
% standard deviation, and runtime mean for the current azimuth angle.
% error_Ment_avg is used as the RMSE-related average error in the simulation
% figures, while sigma_error corresponds to the STD metric in Eq. (73).
S_meas_avg = sum(S_meas_Ment,2)/Ment_Num;
error_Ment_avg = sum(error_Ment)/Ment_Num;
sigma_error = sqrt(sum((error_Ment(1,:)-error_Ment_avg).^2,2)/(Ment_Num-1));
runtime_Ment_avg = sum(runtime_Ment)/Ment_Num;

%% ====================== Outputs for the current azimuth angle ======================
% Store the outputs for the current azimuth angle. These quantities are used
% to compute and visualize the evaluation metrics in Section 4 of the paper.
%
% S_meas_seita(:,seita_iter)     : averaged estimated target position over
%                                 repeated Monte Carlo trials.
% error_avg_seita(seita_iter)    : average localization error over Monte Carlo
%                                 trials, used for RMSE-related analysis in
%                                 Eq. (68).
% sigma_error_seita(seita_iter)  : standard deviation of localization error,
%                                 corresponding to the STD metric in Eq. (73).
% LinHPS_CRLB(seita_iter)        : Cramér-Rao lower bound of the current
%                                 target position, corresponding to Eq. (69).
% runtime_avg_seita(seita_iter)  : average runtime of the online LinHPS
%                                 localization procedure.
S_meas_seita(:,seita_iter) = S_meas_avg;
error_avg_seita(seita_iter) = error_Ment_avg;
sigma_error_seita(seita_iter) = sigma_error;
LinHPS_CRLB(seita_iter) = CRLB(Anchor,MS);
runtime_avg_seita(seita_iter) = runtime_Ment_avg;
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
