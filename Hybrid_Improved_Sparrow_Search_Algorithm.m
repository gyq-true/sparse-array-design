%% 大规模平面稀疏阵列设计：三步走集成算法 (32x32 初始规模)
% 第一步：凸优化 (L1重加权) 优化阵元数目与初步筛选
% 第二步：ISSA 优化阵元位置 (非等参微调)
% 第三步：凸优化 优化最终布局的阵元幅度
%------------------------------------------------
%最终阵元数231，最终PSLL：-22.28dB。工作空间保存为Hybrid_Sparrow_search.mat
clear; clc;

%% --- 参数设置 ---
N_x = 32; N_y = 32;         
N_candidate = N_x * N_y;    % 初始候选位置 1024
lambda = 1;                 
d_init = 0.5 * lambda;      
d_min = 0.35 * lambda;       
L_x = (N_x-1)*d_init; 
L_y = (N_y-1)*d_init; 
[X_grid, Y_grid] = meshgrid(0:d_init:L_x, 0:d_init:L_y);
X_init = X_grid(:);
Y_init = Y_grid(:);

% UV 采样点（高密度采样以保证 -30dB 约束有效）
k = 2 * pi / lambda;
[U, V] = meshgrid(linspace(-1,1,60));
UV_dist = sqrt(U.^2 + V.^2);
UV_dist=UV_dist(:);
u1 = 0.12;               % 主瓣半径
mask_sll = (UV_dist >= u1);
main_mask=(UV_dist <= u1);
u_s = U(mask_sll); v_s = V(mask_sll);%旁瓣区域
A_init = exp(1j * k* (X_init * u_s' + Y_init * v_s'));%旁瓣

%高斯理想响应
sigma = 0.05;           % 控制理想主瓣宽度的参数
F_des = exp(-(UV_dist.^2) / (2 * sigma^2));
F_des_sll = F_des(mask_sll); % 理想响应在副瓣区的分布（接近0）
u_all = U(:); v_all = V(:);
F_des_all = F_des(:);
A_all = exp(1j * k * (X_init * u_all' + Y_init * v_all'));

%% ================== 第一步：凸优化优化阵元数目 ==================
% 迭代重加权 L1 最小化 (Reweighted L1)
w = ones(N_candidate, 1);
weights_l1 = ones(N_candidate, 1);
SLL_target_dB = -20;    % 【设定】旁瓣目标 -30dB
SLL_limit = 10^(SLL_target_dB/20);
w_cvx_opt = ones(N_candidate, 1) / N_candidate; % 保底解
for iter = 1:5
    cvx_begin quiet
        variable w_cvx(N_candidate)
        minimize( norm(weights_l1 .* w_cvx, 1) )
        subject to
            sum( w_cvx ) == 1; % 修正：主瓣峰值归一化
            abs(A_init' * w_cvx) <= SLL_limit; % 注意是 A_init' 
            w_cvx >= 0; w_cvx <= 1; % 物理界限约束

    cvx_end
    
    % 防毒墙：确保只有求解成功才更新权重
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        w_cvx_opt = w_cvx; 
        weights_l1 = 1 ./ (abs(w_cvx_opt) + 0.01);
    else
        fprintf('警告: 凸优化在第 %d 次迭代无解 (可能 -30dB 过于严苛)，停止重加权。\n', iter);
        break;
    end
end
w_cvx = w_cvx_opt; % 将最后一次成功的解赋给 w_cvx 供后续筛选

% 筛选活跃阵元 (保留激励较大的阵元)
active_idx = find(w_cvx > 1e-8);
N_active = length(active_idx);
X_active = X_init(active_idx);
Y_active = Y_init(active_idx);
% 提取保留阵元对应的流形矩阵行
A_init_active = A_init(active_idx, :); 
A_all_active  = A_all(active_idx, :);
fprintf('初步筛选完成，保留阵元数: %d\n', N_active);

%% Tent映射
    % delta_limit: 扰动范围限制
    % Pop_Size: 种群大小 Np
    % dim: 变量维度 (N_active * 2)
    % delta_limit: 扰动范围限制
    delta_limit=lambda/60;
    Np=6;%麻雀的只数
    dim=N_active*2;
    z = zeros(Np, dim);
    % 1. 随机产生初始种子 (避免取 0.5 等不动点)
    z_curr = rand(1, dim); 
    % 2. 混沌迭代生成序列 (公式 15 逻辑简化版)
    for i = 1:Np
        for j = 1:dim
            % Tent 映射结合 Bernoulli 特性的典型实现
            if z_curr(j) < 0.5
                z_curr(j) = 2 * z_curr(j);
            else
                z_curr(j) = 2 * (1 - z_curr(j));
            end
            
            % 引入微小随机扰动避免陷入死循环 (伯努利变换的思想)
            z_curr(j) = mod(z_curr(j) + rand*0.01, 1);
            
            z(i, j) = z_curr(j)+rand/Np^2;
        end
    end
    % 3. 映射到实际搜索空间 [-delta, delta]
    pop_pos = (z - 0.5) * 2 * delta_limit;


%% ================== 第二步：ISSA 优化阵元位置 (严格机制重构) ==================
Max_Iter = 15;
Pop_Size = 6;%麻雀的只数，即阵列布局方案，等于Np

% 1. 初始化种群扰动
fitness_pop = zeros(Np, 1);%每套方案的适应度

% 2. 初始种群预评估 (必须先有适应度，才能在第一次迭代排定身份)
for i = 1:Pop_Size
    pop_pos(i, :) = apply_strict_constraints(pop_pos(i, :), X_active, Y_active, d_min, L_x, L_y, N_active);
    X_curr = X_active + pop_pos(i, 1:N_active)';
    Y_curr = Y_active + pop_pos(i, N_active+1:end)';
    [~, fitness_pop(i)] = solve_convex_amplitudes(X_curr, Y_curr, lambda);
end

% 记录全局最优
[best_fitness, best_idx] = min(fitness_pop);
best_pos_vec = pop_pos(best_idx, :);

% ISSA 角色比例设定，麻雀指代布局方案
PD_num = max(1, round(0.2 * Pop_Size)); % 探险者数量
SD_num = max(1, round(0.2 * Pop_Size)); % 警戒者数量
ST = 0.8; % 安全阈值 (Safety Threshold)

% 3. ISSA 核心主循环
for t = 1:Max_Iter
fprintf('Iteration:%d,正在通过 ISSA 优化非等参位置 \n',Max_Iter);
    % --- [核心机制 1：基于 PSLL 排序分配身份] ---
    [sort_fit, sort_idx] = sort(fitness_pop, 'ascend'); % 升序，越小越好
    pop_pos_sort = pop_pos(sort_idx, :); % 按性能排好序的种群
    
    new_pop_pos = zeros(Pop_Size, dim);   % 更新后的位置
    best_pos_curr = pop_pos_sort(1, :);   % 当前代最佳位置 (探险者领头羊)
    worst_pos_curr = pop_pos_sort(end, :);% 当前代最差位置
    best_f_curr = sort_fit(1);   %最佳的适应度
    worst_f_curr = sort_fit(end);   %最差的适应度
    upper_bound = max(pop_pos_sort);
    lower_bound = min(pop_pos_sort);

% --- [角色 1：探险者 Discoverers (排名前 20%)] ---
for i = 1:PD_num
    R2 = rand; % 预警值
    
    % 1. 计算常规更新位置 (正向位置)
    if R2 < ST 
        % 环境安全：在当前最优解附近大范围广度搜索
        pos_normal = pop_pos_sort(i, :) .* exp(-t / (rand * Max_Iter));
    else
        % 发现危险：随机游走，试图跳出局部最优
        pos_normal = pop_pos_sort(i, :) + randn(1, dim);
    end
    
    % 2. 计算动态反向位置 (Dynamic OBL)
    % 生成关于当前有效搜索空间中心对称的坐标
    pos_reverse = rand(1, dim) .* (upper_bound + lower_bound) - pos_normal; 
    
    % 3. 对两个位置分别进行物理约束修正
    pos_normal = apply_strict_constraints(pos_normal, X_active, Y_active, d_min, L_x, L_y, N_active);
    pos_reverse = apply_strict_constraints(pos_reverse, X_active, Y_active, d_min, L_x, L_y, N_active);
    
    % 4. 还原为实际天线坐标
    X_n = X_active + pos_normal(1:N_active)';
    Y_n = Y_active + pos_normal(N_active+1:end)';
    
    X_r = X_active + pos_reverse(1:N_active)';
    Y_r = Y_active + pos_reverse(N_active+1:end)';
    
    % 5. 当场调用凸优化进行评估 (择优录取)
    [~, fit_normal]  = solve_convex_amplitudes(X_n, Y_n, lambda);
    [~, fit_reverse] = solve_convex_amplitudes(X_r, Y_r, lambda);
    
    % 6. 贪婪保留：谁的副瓣(PSLL)更低，就选谁作为探险者的新位置
    if fit_reverse < fit_normal
        new_pop_pos(i, :) = pos_reverse;
        % 可选：将 fit_reverse 存入一个数组，避免主循环末尾重复计算
    else
        new_pop_pos(i, :) = pos_normal;
        % 可选：存入 fit_normal
    end
end

    % --- [角色 2：跟随者 Followers (排名靠后的 80%)] ---
for i = (PD_num + 1):Pop_Size
    if i > Pop_Size / 2
        % 处于极度饥饿状态 (最差的个体)：飞往其他地方随机搜索
        % 对应公式上子式：Q_r * exp(...)
        new_pop_pos(i, :) = randn(1, dim) .* exp((worst_pos_curr - pop_pos_sort(i, :)) / i^2);
    else
        % 正常的跟随者：向当前的最佳探险者靠拢
        % 对应公式下子式：d_p + |d_i - d_p| * A^+ * L_1
        % 1. 生成 1 x dim 的随机 1 或 -1 矩阵
        A = floor(rand(1, dim) * 2) * 2 - 1; 
        % 2. 计算 A 的伪逆矩阵 (A^+)
        % A * A' 是向量内积，结果是标量 dim，所以伪逆是一个 dim x 1 的列向量
        A_plus = A' / (A * A'); 
        % 3. 生成 1 x dim 的全 1 矩阵 L1
        L1 = ones(1, dim); 
        % 4. 【关键修正】执行严格的矩阵连乘 (*)
        % 计算结果：abs(...) * A_plus 会算出一个标量，再乘以 L1 变成全维度相同步长的行向量
        new_pop_pos(i, :) = best_pos_curr + abs(pop_pos_sort(i, :) - best_pos_curr) * A_plus * L1;
    end
end

    % --- [角色 3：警戒者 Scouters (随机抽取 20%)] ---
    % 警戒者是在探险者和跟随者更新完基础位置后，随机挑选几只产生避险行为
    scout_idx = randperm(Pop_Size, SD_num);
    for j = 1:SD_num
        idx = scout_idx(j);
        if sort_fit(idx) > best_f_curr
            % 自身不是最优：向全局最优位置靠拢
            new_pop_pos(idx, :) = best_pos_curr + randn(1, dim) .* abs(new_pop_pos(idx, :) - best_pos_curr);
        else
            % 自身已经是全局最优 (意识到高处不胜寒的危险)：在自身附近微小游走
            new_pop_pos(idx, :) = new_pop_pos(idx, :) + (randn(1, dim) * 2 - 1) .* abs(new_pop_pos(idx, :) - worst_pos_curr) / (sort_fit(idx) - worst_f_curr + 1e-8);
        end
    end
% =====================================================================
    % --- [角色 4：自适应 T 分布变异扰动 (结合收敛概率机制)] ---
    % 核心思路：按一定概率对“当前最优解”进行变异，并将变异后的图纸交给“当前最差个体”去试错
    % =====================================================================
    
    % 1. 计算后期收敛概率 (动态变异概率)
    % 采用非线性递减概率：前期概率大（鼓励变异探索），后期概率小（收敛于最优解）
    p_mutation = exp(-5 * t / Max_Iter); 
    
    if rand < p_mutation
        % 2. 核心数学：生成 T 分布扰动向量
        % 自由度 (Degrees of freedom) 严格设定为当前迭代次数 t
        % trnd 是 MATLAB 统计工具箱中生成 T 分布随机数的自带函数
        t_step = trnd(t, 1, dim); 
        
        % 3. 变异公式：在新位置中注入 T 分布扰动
        % 注意：扰动是乘性的，且直接作用于 best_pos_curr
        mutated_pos = best_pos_curr + best_pos_curr .* t_step;
        
        % 4. 零成本算力替换（极其关键的工程技巧）
        % 为了不增加额外的 CP 调用次数，我们将这份“变异图纸”强行塞给本代最差的麻雀
        % new_pop_pos 的最后一行 (end) 恰好对应本代排位倒数第一的麻雀
        new_pop_pos(end, :) = mutated_pos; 
        
        fprintf('  -> 触发 T 分布变异: 最差个体已被最优解的变异分支覆盖。\n');
    end
    
    % --- [核心机制 2：物理约束检查与 CP 重新评估] ---
    for i = 1:Pop_Size
        % 1. 执行严格的边界与最小距离约束
        new_pop_pos(i, :) = apply_strict_constraints(new_pop_pos(i, :), X_active, Y_active, d_min, L_x, L_y, N_active);
        % 2. 获取实际天线物理坐标
        X_curr = X_active + new_pop_pos(i, 1:N_active)';
        Y_curr = Y_active + new_pop_pos(i, N_active+1:end)';
        
        % 3. CP 求解该非等参布局下的最优幅度与 PSLL
        [~, current_psll] = solve_convex_amplitudes(X_curr, Y_curr, lambda);
        
        % 4. 贪婪更新 (仅当新位置比老位置好时才接受，防止物理约束导致位置变差)
        if current_psll < fitness_pop(sort_idx(i))
            pop_pos(sort_idx(i), :) = new_pop_pos(i, :);
            fitness_pop(sort_idx(i)) = current_psll;
        end
    end
    
    % --- [记录历史全局最优] ---
    [current_best_fit, min_idx] = min(fitness_pop);
    if current_best_fit < best_fitness
        best_fitness = current_best_fit;
        best_pos_vec = pop_pos(min_idx, :);
    end
    
    fprintf('ISSA 迭代 %d: 当前最佳 PSLL = %.2f dB\n', t, best_fitness);
end

% 提取并固定最终位置
X_final_pos = X_active + best_pos_vec(1:N_active)';
Y_final_pos = Y_active + best_pos_vec(N_active+1:end)';

%% ================== 第三步：凸优化优化阵元幅度 ==================
fprintf('步骤 3: 正在对最终布局进行高精度幅度加权优化...\n');
% 调用高精度版本，步长从 4 降到 1 或 0.5
[final_weights, final_psll] = solve_convex_amplitudes_fine(X_final_pos, Y_final_pos, lambda);

%% 结果展示
fprintf('\n--- 最终设计结果 ---\n');
fprintf('最终阵元数: %d\n', N_active);
fprintf('最终 PSLL: %.2f dB\n', final_psll);

figure;
subplot(1,2,1);
scatter(X_final_pos, Y_final_pos, 25, final_weights, 'filled');
title('最终非等参稀疏布局 ');
xlabel('x (\lambda)'); ylabel('y (\lambda)'); axis equal; colorbar;

subplot(1,2,2);
% 绘制方向图切面 (Phi = 0)
theta_plot = -90:0.5:90;
u_plot = sin(deg2rad(theta_plot));
AF = abs(exp(1j * (2*pi/lambda) * (X_final_pos * u_plot))' * final_weights);
plot(theta_plot, 20*log10(AF/max(AF)), 'LineWidth', 1.5);
title('方向图切面 (\phi=0^\circ)');
xlabel('\theta (deg)'); ylabel('Normalized Pattern (dB)');
grid on; ylim([-60 0]);
%% %%  追加计算：半功率主瓣宽度 (HPBW) 与 增益 ---
fprintf('\n--- 天线核心指标分析 ---\n');

% ==========================================
% 1. 计算 HPBW (在 phi = 0 度切面上)
% ==========================================
% 为了保证精度，在主瓣附近进行极高密度采样 (步长 0.001 度)
theta_hpbw = -10:0.001:10; 
u_hpbw = sin(deg2rad(theta_hpbw));
AF_hpbw = abs(exp(1j * (2*pi/lambda) * (X_final_pos * u_hpbw))' * final_weights);
AF_hpbw_dB = 20 * log10(AF_hpbw / max(AF_hpbw));

% 寻找峰值位置 (理论上在 0 度)
[~, peak_idx] = max(AF_hpbw_dB);

% 在峰值左侧寻找最接近 -3dB 的点
[~, left_idx] = min(abs(AF_hpbw_dB(1:peak_idx) - (-3)));
theta_left = theta_hpbw(left_idx);

% 在峰值右侧寻找最接近 -3dB 的点
[~, right_idx] = min(abs(AF_hpbw_dB(peak_idx:end) - (-3)));
right_idx = right_idx + peak_idx - 1;
theta_right = theta_hpbw(right_idx);

% 最终半功率主瓣宽度
HPBW = theta_right - theta_left;
fprintf('半功率主瓣宽度 (HPBW): %.3f° (从 %.3f° 到 %.3f°)\n', HPBW, theta_left, theta_right);

% ==========================================
% 2. 计算主瓣增益与阵列方向性 (Directivity)
% ==========================================


% 【学术常用】通过空间球面积分计算真正的阵列方向性 (Directivity)
% 设置球面积分的步长 (1度足够保证精度)
d_theta = 1; 
d_phi = 1;
theta_int = 0:d_theta:90; 
phi_int = 0:d_phi:359;
[TH, PH] = meshgrid(deg2rad(theta_int), deg2rad(phi_int));

% 将球坐标转化为 UV 坐标
U_int = sin(TH) .* cos(PH);
V_int = sin(TH) .* sin(PH);

% 计算整个半球面所有积分点上的 Array Factor (AF) 的平方
AF_int = abs(exp(1j * (2*pi/lambda) * (X_final_pos * U_int(:)' + Y_final_pos * V_int(:)'))' * final_weights).^2;
AF_int = reshape(AF_int, length(phi_int), length(theta_int));

% 计算半球面辐射总功率 Prad_half (数值积分)
Prad_half = sum(sum(AF_int .* sin(TH) * deg2rad(d_theta) * deg2rad(d_phi)));

% 全空间辐射总功率 (由于平面阵列的波束前后对称，全空间总功率 = 半空间的 2 倍)
Prad_full = 2 * Prad_half;

% 峰值辐射强度 U_max
U_max = AF_peak_raw^2;

% 计算方向性 Directivity = 4*pi*U_max / Prad_full
Directivity_linear = 4 * pi * U_max / Prad_full;
Directivity_dBi = 10 * log10(Directivity_linear);

fprintf('最终阵列方向性 (Directivity): %.2f dBi\n', Directivity_dBi);


%% ================== 局部函数区 ==================

function [w, psll] = solve_convex_amplitudes(X, Y, lambda)
    N = length(X);
    k = 2 * pi / lambda;
    theta_limit = 4;
    [ts, ps] = meshgrid(theta_limit+1:4:90, 0:45:315);
    us = sin(deg2rad(ts(:))) .* cos(deg2rad(ps(:)));
    vs = sin(deg2rad(ts(:))) .* sin(deg2rad(ps(:)));
    As = exp(1j * k * (X .* us' + Y .* vs'));
    As_H = As'; % 在 CVX 外部提前做共轭转置
    try
        cvx_begin quiet
            variable wv(N)
            variable g
            minimize( g )
            subject to
                sum(wv) == 1;
                abs(As_H * wv) <= g;
                wv >= 0; wv <= 1;
        cvx_end
        psll = 20 * log10(double(g));
        w = double(wv);
    catch
        w = ones(N, 1)/N; psll = -10;
    end
end

function nvec = apply_strict_constraints(pvec, Xa, Ya, dm, Lx, Ly, Na)
%pop_pos(i, :) = apply_strict_constraints(pop_pos(i, :), X_active, Y_active, d_min, L_x, L_y, N_active);
    dx = pvec(1:Na)'; 
    dy = pvec(Na+1:end)';
    X = max(0, min(Lx, Xa + dx)); % 第一次边界截断
    Y = max(0, min(Ly, Ya + dy));
    % 设置最大排斥迭代次数（防止死循环，通常3-5次足以解开大部分密集碰撞）
    max_repulsion_iters = 5; 
    for iter = 1:max_repulsion_iters
        collision_found = false; % 碰撞标记
        for i = 1:Na
            for j = i+1:Na
                dist = sqrt((X(i)-X(j))^2 + (Y(i)-Y(j))^2);
                if dist < dm
                    collision_found = true;
                    % 1. 计算排斥方向向量 (从 i 指向 j)
                    if dist < 1e-6
                        % 极端情况：两点完全重合，随机给一个方向弹开
                        angle = rand * 2 * pi;
                        dir_x = cos(angle);
                        dir_y = sin(angle);
                    else
                        % 归一化方向向量
                        dir_x = (X(j) - X(i)) / dist;
                        dir_y = (Y(j) - Y(i)) / dist;
                    end
                    % 2. 计算需要推开的重叠距离 (加一点微小裕度 1e-4 确保彻底分开)
                    overlap = dm - dist + 1e-8; 
                    half_overlap = overlap / 2; % 两人各退一步，分摊位移
                    % 3. 沿连线方向反向推开
                    X(i) = X(i) - dir_x * half_overlap;
                    Y(i) = Y(i) - dir_y * half_overlap;
                    X(j) = X(j) + dir_x * half_overlap;
                    Y(j) = Y(j) + dir_y * half_overlap;
                    % 4. 立即进行边界安全截断 (防止被挤出孔径)
                    X(i) = max(0, min(Lx, X(i)));
                    Y(i) = max(0, min(Ly, Y(i)));
                    X(j) = max(0, min(Lx, X(j)));
                    Y(j) = max(0, min(Ly, Y(j)));
                end
            end
        end
        % 如果这一轮历遍没有发现任何碰撞，说明全部安全，提前结束迭代
        if ~collision_found
            break;
        end
    end
    % 最后将绝对坐标重新转化为偏移量返回
    nvec = [(X - Xa)', (Y - Ya)'];
end
function [w, psll] = solve_convex_amplitudes_fine(X, Y, lambda)
    N = length(X);
    k = 2 * pi / lambda;
    theta_limit = 4;
    
    % 【关键改变】：极大地加密采样点！
    % theta 步长设为 1度，phi 步长设为 10度 (原来是 4度 和 45度)
    [ts, ps] = meshgrid(theta_limit+1 : 1 : 90, 0 : 10 : 359); 
    
    us = sin(deg2rad(ts(:))) .* cos(deg2rad(ps(:)));
    vs = sin(deg2rad(ts(:))) .* sin(deg2rad(ps(:)));
    As = exp(1j * k * (X .* us' + Y .* vs'));
    As_H = As'; % 在 CVX 外部提前做共轭转置
    try
        cvx_begin quiet
            variable wv(N)
            variable g
            minimize( g )
            subject to
                sum(wv) == 1;
                abs(As_H * wv) <= g;
                wv >= 0; wv <= 1;
        cvx_end
        psll = 20 * log10(double(g));
        w = double(wv);
    catch
        w = ones(N, 1)/N; psll = -10;
    end
end

