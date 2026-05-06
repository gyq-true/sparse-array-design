%% ============================================================
%  基于 COADBO 的平面稀疏阵列综合（复现论文）
%  计算适应度的函数里方向性稀疏可以写到目标里吗？
% ============================================================
% 比较好的一次运行结果保存为method2.mat，psll=-20.982dB，密度锥化阵元分布后为512个阵元。主瓣增益为27.66dBi
% 主瓣宽度为4.634度
% ============================================================
clear; clc;
 
%% 参数设置（与论文一致）
N_x = 32; N_y = 32;   % 论文示例规模
lambda = 1;
d = 0.5 * lambda;
k = 2*pi/lambda;

L_x = (N_x-1)*d;
L_y = (N_y-1)*d;

[X_grid, Y_grid] = meshgrid(0:d:L_x, 0:d:L_y);
X_init = X_grid(:);
Y_init = Y_grid(:);
N = length(X_init);


PSLL_d = -25;           % 目标峰值旁瓣电平 (dB)
DIR_d = 27.32;           % 目标方向性系数 (dBi)
dim = 2 * N;

[U, V] = meshgrid(linspace(-1,1,60));
UV_dist = sqrt(U.^2 + V.^2);
U_grid=U(:);V_grid=V(:);
UV_dist=UV_dist(:);
u1 = 0.12;               % 主瓣半径
mask_sll = (UV_dist >= u1);
main_mask=(UV_dist <= u1);%主瓣区域
u_s = U(mask_sll); v_s = V(mask_sll);%旁瓣区域
%% COADBO 参数
Pop_Size = 10;
Max_Iter = 30;

%计算上下界
L = 15.5; H = 15.5;
lb = zeros(1, dim);
ub = [repmat(L, 1, N), repmat(H, 1, N)];
%% ================== 混沌初始化 + OBL ==================
% --- Fuch混沌映射 ---

pop_pos = zeros(Pop_Size, dim);
z = rand(1, dim);
for i = 1:Pop_Size
    z = cos(1 ./ (z.^2 + 1e-6)); % Fuch 映射简化形式
    pop_pos(i,:) = lb + (ub - lb) .* abs(mod(z, 1));
end

%% ================== 初始适应度 ==================
fitness = zeros(Pop_Size,1);
fitness_obl = zeros(Pop_Size, 1);

for i = 1:Pop_Size
    X =  pop_pos(i,1:N)';
    Y =  pop_pos(i,N+1:end)';
    [fitness(i), ~] = solve_convex_amplitudes(X, Y, lambda,PSLL_d,DIR_d);
end

%% ================== COADBO 主循环 ==================
for t = 1:Max_Iter
    for i = 1:Pop_Size
        pos_vec = pop_pos(i, :);
        qx =  pos_vec(1:N)';
        qy =  pos_vec(N+1:end)';

        try
            [current_fitness, current_psll, current_dir] = solve_convex_amplitudes(qx, qy, lambda,PSLL_d,DIR_d);
        catch
            current_psll = 0; current_dir = 0; % 容错
        end
        
        if current_fitness < fitness(i)
            pop_pos(i,:) = [qx',qy'];
            fitness(i) = current_fitness;
        end
        K = rand;
        % 必须使用全局的上下界 lb 和 ub (维度为 1×dim)
        pop_pos_obl = lb + ub - K *pop_pos(i,:); 
        % 边界越界保护
        pop_pos_obl = max(lb, min(ub, pop_pos_obl));
        % 计算逆解种群的适应度
        pos_vec_i = pop_pos_obl;
        qxi = pos_vec_i(1:N);
        qyi = pos_vec_i(N+1:end);
        
        try
            [fitness_obl(i), ~, ~, ~] = solve_convex_amplitudes(qxi', qyi', lambda, PSLL_d, DIR_d);
        catch
            fitness_obl(i) = 1e9;
        end
        if fitness(i) < fitness_obl(i)
            pop_pos(i,:) = [qx',qy'];%原适应度更好
        else%逆解适应度更好
            pop_pos(i,:) = [qxi,qyi];
            fitness(i) = fitness_obl(i);
        end
    end

   
         %% 排序与最优位置更新
    [fitness, sorted_idx] = sort(fitness);
    pop_pos = pop_pos(sorted_idx, :);
    best_pos = pop_pos(1, :);
    worst_pos = pop_pos(end, :);
    
    % 分配蜣螂角色 (比例参考标准DBO)
    n_ball = round(Pop_Size * 0.2);     % 滚球
    n_breed = round(Pop_Size * 0.2);    % 繁殖
    n_forage = round(Pop_Size * 0.35);  % 觅食
    n_steal = Pop_Size - n_ball - n_breed - n_forage; % 偷窃
    
    new_pop = pop_pos;
    
    % (1) 滚球阶段：引入 OOA 全球探索策略 (公式 6) 
    eta = randi([1, 2]); % 随机数 η ∈ {1, 2} [cite: 537]
    for i = 1:n_ball
        SF = pop_pos(randi(Pop_Size), :); % 随机选择食物
        new_pop(i, :) = pop_pos(i, :) + rand(1, dim) .* (SF - eta * pop_pos(i, :));
    end
    
    % (2) 繁殖阶段：标准 DBO 边界收缩 [cite: 564]
    R = 1 - t/Max_Iter;
    % 根据文献中觅食阶段的公式：基于最优解的下界和上界
    lb_b = max(best_pos * (1 - R),lb); 
    ub_b = max(best_pos * (1 + R),ub);
    for i = n_ball+1 : n_ball+n_breed
        new_pop(i, :) = best_pos + rand(1, dim) .* (rand(1, dim).*(lb_b - ub_b));
    end
    
    % (3) 觅食阶段：自适应步长策略 (公式 8) [cite: 550, 565]
    alpha_0 = cos(pi/3 * (1 + t/Max_Iter)); % 自适应步长因子 [cite: 549]
    for i = n_ball+n_breed+1 : n_ball+n_breed+n_forage
        C1 = rand(1, dim); C2 = rand(1, dim);
        new_pop(i, :) = alpha_0 * pop_pos(i, :) + C1.*(pop_pos(i,:) - lb_b) + C2.*(pop_pos(i,:) - ub_b);
    end
    
    

    % (4) 偷窃阶段：自适应 t-分布扰动 (公式 9)
    r_deg = exp(1 + (t/Max_Iter)^2); % 自由度 r
    for i = Pop_Size-n_steal+1 : Pop_Size
        % 1. 生成标准正态分布随机数 X ~ N(0, 1)
        norm_x = randn(1, dim);
        % 2. 模拟自由度为 r_deg 的卡方分布随机数
        % chisq_y 等价于 r_deg 个独立标准正态分布的平方和
        chisq_y = sum(randn(floor(r_deg), dim).^1, 1); % 这里通过近似或直接循环求和实
        % 使用更标准、直接的卡方分布模拟方法 (Gamma 分布特例)
        % 这里我们用标准正态分布的平方和来近似或计算卡方分布
        Z_norm = randn(round(r_deg), dim);
        chisq_r = sum(Z_norm.^2, 1);
        
        % 3. 计算 t-分布随机数
        t_dist = norm_x ./ sqrt(chisq_r / r_deg); 
        
        new_pop(i, :) = best_pos + best_pos .* t_dist;
    end

        %% --- 边界控制 ---
        pop_pos = max(lb, min(ub, pop_pos));

    best_qx = pop_pos(1, 1:N);
    best_qy = pop_pos(1, N+1:end);
    
    [best_fitness,  best_psll, best_dir] = solve_convex_amplitudes(best_qx', best_qy', lambda, PSLL_d, DIR_d);
    
    fprintf('Iter %d: Best Obj = %.4f | PSLL = %.2f dB | Directivity = %.2f dBi\n', ...
            t, best_fitness, best_psll, best_dir);
end
%% 密度锥化法得到稀疏阵列
% ================== 密度锥化与阵元削减 (Post-Processing) ==================

% 2. 计算各阵元距离中心的径向距离 (Normalized Radius)
% 假设中心在 (L/2, H/2)
center_x = L_x/2; center_y = L_y/2;
rho = sqrt((best_qx - center_x).^2 + (best_qy - center_y).^2);
rho_max = max(rho);

% 3. 密度锥化函数分配 (此处采用典型的 Cosine-Squared 锥化作为示例)
% 距离中心越近，权值越高；距离越远，权值越低
weights_tapered = (cos(pi/2 * rho / rho_max)).^2 + 0.1; 

% 4. 阵元削减 (Thinning): 削减 50%
% 按照权值大小排序，只保留权值最大的前 50% 的阵元
num_keep = round(N * 0.5); 
[~, sort_idx] = sort(weights_tapered, 'descend');
keep_idx = sort_idx(1:num_keep);

% 得到削减后的稀疏布局
sparse_qx = best_qx(keep_idx);
sparse_qy = best_qy(keep_idx);
sparse_weights = weights_tapered(keep_idx); % 或者使用 1 (均匀权值)

% 5. 重新计算削减后的性能指标 (使用您之前校准的积分法)
[final_fitness,  final_psll, final_dir] = solve_convex_amplitudes_fine(sparse_qx', sparse_qy', lambda, PSLL_d, DIR_d,sparse_weights);

% ================== 结果输出与对比 ==================
fprintf('\n================ 稀疏削减结果 =================\n');
fprintf('初始阵元数: %d\n', N);
fprintf('削减后阵元数: %d (稀疏率: 50%%)\n', num_keep);
fprintf('最终 PSLL: %.2f dB\n', final_psll);
fprintf('最终方向性系数: %.2f dBi\n', final_dir);
fprintf('最终适应度 Obj: %.4f\n', final_fitness);

%% %% ============================================================
%  基于削减后稀疏阵列的性能分析与可视化
% ============================================================

fprintf('\n--- 最终设计结果（稀疏削减后）---\n');
fprintf('最终稀疏阵元数: %d\n', num_keep);
fprintf('最终 PSLL: %.2f dB\n', final_psll);

% 1. 可视化布局与方向图切面
figure('Name', '稀疏阵列综合结果', 'Color', 'white');

% % 子图 1：阵列几何分布图
% subplot(1, 2, 1);
% scatter(sparse_qx, sparse_qy, 35, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
% title('削减后的稀疏平面布局 (50% 稀疏率)');
% xlabel('x (\lambda)'); ylabel('y (\lambda)'); 
% axis equal; grid on; box on;
% 子图 1：阵列几何分布图
subplot(1, 2, 1);
% 关键修改：将固定颜色替换为 sparse_weights，让颜色反映幅度
scatter(sparse_qx, sparse_qy, 35, sparse_weights, 'filled'); 

% 添加色带并设置标签
colormap(gca, 'jet'); % 推荐使用 jet 或 parula 色系
cb = colorbar; 
cb.Label.String = 'Excitation Amplitude (幅度)'; % 给 colorbar 加上说明
cb.Label.FontSize = 10;

title(sprintf('削减后的稀疏平面布局 (保留阵元: %d)', num_keep));
xlabel('x (\lambda)'); ylabel('y (\lambda)'); 
axis equal; grid on; box on;

% 子图 2：方向图切面 (\phi=0^\circ)
subplot(1, 2, 2);
theta_plot = -90:0.5:90;
u_plot = sin(deg2rad(theta_plot));
% 计算削减后的方向图切面
u_plot_row = u_plot(:)'; 
qx_col = sparse_qx(:);
weights = sparse_weights(:);
AF = abs(weights.' * exp(1j * k * (qx_col * u_plot_row)));
plot(theta_plot, 20*log10(AF/max(AF)), 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
title('方向图切面 (\phi=0^\circ)');
xlabel('\theta (deg)'); ylabel('Normalized Pattern (dB)');
grid on; ylim([-50 0]); xlim([-90 90]);

% 2. 追加计算：半功率主瓣宽度 (HPBW)
fprintf('\n--- 天线核心指标分析 ---\n');

theta_hpbw = -10:0.001:10; 
u_hpbw_row = sin(deg2rad(theta_hpbw(:)'));
AF_hpbw = abs(weights.' * exp(1j * k * (qx_col * u_hpbw_row)));
AF_hpbw_dB = 20 * log10(AF_hpbw / max(AF_hpbw));

% 寻找峰值位置
[AF_peak_raw, peak_idx] = max(AF_hpbw);
[~, left_idx] = min(abs(AF_hpbw_dB(1:peak_idx) - (-3)));
theta_left = theta_hpbw(left_idx);

[~, right_idx] = min(abs(AF_hpbw_dB(peak_idx:end) - (-3)));
right_idx = right_idx + peak_idx - 1;
theta_right = theta_hpbw(right_idx);

HPBW = theta_right - theta_left;
fprintf('半功率主瓣宽度 (HPBW): %.3f° (从 %.3f° 到 %.3f°)\n', HPBW, theta_left, theta_right);

% 3. 追加计算：真实阵列方向性 (Directivity)
d_theta = 1; 
d_phi = 1;
theta_int = 0:d_theta:90; 
phi_int = 0:d_phi:359;
[TH, PH] = meshgrid(deg2rad(theta_int), deg2rad(phi_int));

U_int = sin(TH) .* cos(PH);
V_int = sin(TH) .* sin(PH);
qx_col = sparse_qx(:);
qy_col = sparse_qy(:);
U_vec = U_int(:)'; 
V_vec = V_int(:)';
AF_vec = abs(sparse_weights(:).' * exp(1j * k * (qx_col * U_vec + qy_col * V_vec)));
% 计算半球面所有积分点上的 Array Factor 的平方
AF_int = reshape(AF_vec.^2, length(phi_int), length(theta_int));

% 计算峰值辐射强度 U_max (注意：已是 AF_int 中的平方最大值)
[U_max, ~] = max(AF_int(:)); 

% 数值积分计算半球面辐射功率
Prad_half = sum(sum(AF_int .* sin(TH) * deg2rad(d_theta) * deg2rad(d_phi)));
Prad_full = 2 * Prad_half;

% 计算方向性
Directivity_linear = 4 * pi * U_max / Prad_full;
Directivity_dBi = 10 * log10(Directivity_linear);

fprintf('最终阵列方向性 (Directivity): %.2f dBi\n', Directivity_dBi);
%% % --- 追加：3D 方向图 u-v 平面可视化（找出隐藏的高旁瓣） ---
figure('Name', 'u-v 平面全局方向图', 'Color', 'white');

% 生成高密度的 u-v 网格
u_range = linspace(-1, 1, 300);
v_range = linspace(-1, 1, 300);
[uu, vv] = meshgrid(u_range, v_range);

% 提取物理可见区域 (u^2 + v^2 <= 1)
visible_region = (uu.^2 + vv.^2 <= 1);

% 重新准备列向量
qx_col = sparse_qx(:);
qy_col = sparse_qy(:);
weights = sparse_weights(:);

% 矩阵计算全局 3D 阵列因子 (只算可见区能极大提速)
AF_uv = zeros(size(uu));
uu_vis = uu(visible_region)';
vv_vis = vv(visible_region)';
AF_vis = abs(weights.' * exp(1j * k * (qx_col * uu_vis + qy_col * vv_vis)));

% 归一化并转 dB
AF_uv(visible_region) = 20 * log10(AF_vis / max(AF_vis));
AF_uv(~visible_region) = NaN; % 可见区外置为空白

% 绘制 2D 伪彩图
contourf(uu, vv, AF_uv, -40:2:0, 'LineStyle', 'none');
colormap('jet');
c = colorbar;
c.Label.String = 'Normalized Directivity (dB)';
hold on;

% 标出主瓣和 -13dB 旁瓣的等高线以示醒目
contour(uu, vv, AF_uv, [-13.14, -13.14], 'LineColor', 'k', 'LineWidth', 1.5);

title('全空间方向图 (u-v 投影)');
xlabel('u = sin(\theta)cos(\phi)'); 
ylabel('v = sin(\theta)sin(\phi)');
axis equal; axis([-1 1 -1 1]); grid on;


%% 局部函数区
% function [fitness, psll, dir_val] = solve_convex_amplitudes(X, Y, lambda, PSLL_d, DIR_d)
%     X = X(:); 
%     Y = Y(:);
%     k = 2 * pi / lambda;
%     N = length(X);
% 
%     theta_limit = 4;
%     d_theta = 1; d_phi = 4;
%     [ts, ps] = meshgrid(0:d_theta:90, 0:d_phi:359);
%     us = sin(deg2rad(ts(:))) .* cos(deg2rad(ps(:)));
%     vs = sin(deg2rad(ts(:))) .* sin(deg2rad(ps(:)));
%     U_grid=us(:)';
%     V_grid=vs(:)';
%     
%     As = abs(sum(exp(1j * k * (X * U_grid + Y * V_grid)),1));
%     AF_max = max(As);
%     AF_norm = As / AF_max;
% 
%      
%     [M, I] = max(AF_norm(:));
%     AF_norm_psllfind=AF_norm;
%     AF_norm_psllfind(I) = -inf; 
%     [pks, locs] = findpeaks(AF_norm_psllfind, 'MinPeakDistance', 20); 
%     if isempty(pks)
%         max_sll_val = max(AF_norm_psllfind);
%     else
%         max_sll_val = max(pks);
%     end
%     psll = 20 * log10(max_sll_val);
%         
%     % 3. 计算当前阵列的方向性系数 (Directivity)
%     % 计算半球面所有积分点上的 Array Factor 的平方（维度：32760 x 1）
%     AF_int = As.^2;
%     
%     % 求最大辐射强度 AF_peak_raw
%     AF_peak_raw = max(AF_int); 
%     
%     % 数值积分计算半球面辐射功率
%     sin_TH = sin(deg2rad(ts(:)))';
%     Prad_half = sum(AF_int .* sin_TH) * deg2rad(d_theta) * deg2rad(d_phi);
%     
%     % 全空间辐射总功率
%     Prad_full = 2 * Prad_half;
%     
%     % 计算方向性系数
%     Directivity_linear = 4 * pi * AF_peak_raw / Prad_full;
%     dir_val = 10 * log10(Directivity_linear);
%     
%     % 4. 计算目标函数
%     term1 = abs(psll - PSLL_d) / abs(PSLL_d);
%     term2 = abs(dir_val - DIR_d) / DIR_d;
%     
%     fitness = term1 + term2;
% end
% 
% function [fitness, psll, dir_val] = solve_convex_amplitudes_fine(X, Y, lambda, PSLL_d, DIR_d, w)
%     X=X(:);Y=Y(:);
%     w=w(:);
%     N = length(X);
%     k = 2 * pi / lambda;
%     
%     theta_limit = 4;
%     d_theta = 1; d_phi = 4;
%     [ts, ps] = meshgrid(0:d_theta:90, 0:d_phi:359);
%     us = sin(deg2rad(ts(:))) .* cos(deg2rad(ps(:)));
%     vs = sin(deg2rad(ts(:))) .* sin(deg2rad(ps(:)));
%     U_grid=us(:)';
%     V_grid=vs(:)';
%     
%     As = abs(w.'*exp(1j * k * (X * U_grid + Y * V_grid)));
%     AF_max = max(As);
%     AF_norm = As / AF_max;
%     
%       
%     [~, I] = max(AF_norm(:));
%     AF_norm_psllfind=AF_norm;
%     AF_norm_psllfind(I) = -inf; 
%     [pks, locs] = findpeaks(AF_norm_psllfind, 'MinPeakDistance', 20); 
%     if isempty(pks)
%         max_sll_val = max(AF_norm_psllfind);
%     else
%         max_sll_val = max(pks);
%     end
%     psll = 20 * log10(max_sll_val);
%         
%     % 3. 计算当前阵列的方向性系数 (Directivity)
%     % 引入权重 w 计算
%     AF_int = As.^2;
%     
%     % 求最大辐射强度 AF_peak_raw
%     AF_peak_raw = max(AF_int); 
%     
%     % 数值积分计算半球面辐射功率
%     sin_TH = sin(deg2rad(ts(:)))';
%     Prad_half = sum(AF_int .* sin_TH) * deg2rad(4) * deg2rad(45);
%     
%     % 全空间辐射总功率
%     Prad_full = 2 * Prad_half;
%     
%     % 计算方向性系数
%     Directivity_linear = 4 * pi * AF_peak_raw / Prad_full;
%     dir_val = 10 * log10(Directivity_linear);
%     
%     % 4. 计算目标函数
%     term1 = abs(psll - PSLL_d) / abs(PSLL_d);
%     term2 = abs(dir_val - DIR_d) / DIR_d;
%     
%     fitness = term1 + term2;
% end
function [fitness, psll, dir_val] = solve_convex_amplitudes(X, Y, lambda, PSLL_d, DIR_d)
    % 1. 强制转为列向量
    X = X(:); 
    Y = Y(:);
    k = 2 * pi / lambda;
    
    % 2. 建立包含 0 度主瓣的完整网格 
    % 步长：theta 用 1 度保证不漏峰，phi 用 4 度加快运算
    d_theta = 1; d_phi = 4;
    [ts, ps] = meshgrid(0:d_theta:90, 0:d_phi:359);
    us = sin(deg2rad(ts(:))) .* cos(deg2rad(ps(:)));
    vs = sin(deg2rad(ts(:))) .* sin(deg2rad(ps(:)));
    U_grid = us(:)';
    V_grid = vs(:)';
    
    % 3. 矩阵法计算阵列因子 As (1 x M 向量)
    % 等幅激励，相当于权重全是 1
    As = abs(ones(1, length(X)) * exp(1j * k * (X * U_grid + Y * V_grid)));
    
    % 提取真实的主瓣峰值 (此时必然在 theta=0 附近)
    AF_max = max(As);
    AF_norm = As / AF_max;
     
    % 4. 使用逻辑掩码提取 PSLL (安全、快速，不会因为 2D 拍扁产生假峰值)
    theta_limit = 4; % 设定主瓣保护区为 4 度
    is_sidelobe = ts(:)' > theta_limit; % 生成旁瓣区域的布尔索引
    psll = 20 * log10(max(AF_norm(is_sidelobe))); % 直接取旁瓣区域的最大值
        
    % 5. 计算真实方向性系数 (Directivity)
    AF_int = As.^2;
    AF_peak_raw = max(AF_int); 
    
    sin_TH = sin(deg2rad(ts(:)))';
    Prad_half = sum(AF_int .* sin_TH) * deg2rad(d_theta) * deg2rad(d_phi);
    Prad_full = 2 * Prad_half;
    
    Directivity_linear = 4 * pi * AF_peak_raw / Prad_full;
    dir_val = 10 * log10(Directivity_linear);
    
    % 6. 计算目标函数
    term1 = abs(psll - PSLL_d) / abs(PSLL_d);
    term2 = abs(dir_val - DIR_d) / DIR_d;
    fitness = term1 + term2;
end

function [fitness, psll, dir_val] = solve_convex_amplitudes_fine(X, Y, lambda, PSLL_d, DIR_d, w)
    % 1. 强制转为列向量
    X = X(:); Y = Y(:); w = w(:);
    k = 2 * pi / lambda;
    
    % 2. 建立包含 0 度的网格
    d_theta = 1; d_phi = 4;
    [ts, ps] = meshgrid(0:d_theta:90, 0:d_phi:359);
    us = sin(deg2rad(ts(:))) .* cos(deg2rad(ps(:)));
    vs = sin(deg2rad(ts(:))) .* sin(deg2rad(ps(:)));
    U_grid = us(:)';
    V_grid = vs(:)';
    
    % 3. 左乘权重矩阵完成叠加
    As = abs(w.' * exp(1j * k * (X * U_grid + Y * V_grid)));
    
    AF_max = max(As);
    AF_norm = As / AF_max;
    
    % 4. 逻辑掩码提取 PSLL
    theta_limit = 4;
    is_sidelobe = ts(:)' > theta_limit;
    psll = 20 * log10(max(AF_norm(is_sidelobe)));
        
    % 5. 计算真实方向性系数
    AF_int = As.^2;
    AF_peak_raw = max(AF_int); 
    
    sin_TH = sin(deg2rad(ts(:)))';
    Prad_half = sum(AF_int .* sin_TH) * deg2rad(d_theta) * deg2rad(d_phi);
    Prad_full = 2 * Prad_half;
    
    Directivity_linear = 4 * pi * AF_peak_raw / Prad_full;
    dir_val = 10 * log10(Directivity_linear);
    
    % 6. 计算目标函数
    term1 = abs(psll - PSLL_d) / abs(PSLL_d);
    term2 = abs(dir_val - DIR_d) / DIR_d;
    fitness = term1 + term2;
end
