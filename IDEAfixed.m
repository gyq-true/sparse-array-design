%% IDEA Algorithm for Sparse Planar Array Synthesis - Robust Version
% =============================================================
%15轮迭代，阵元数目252，方向维半功率波束宽度4.7043°
% =============================================================
%2026.3.23修改：️
% 1.优化目标：L2,1 混合范数（不好?）/1.5范数(√），结果与IDEA类似，迭代轮次少
% 2.物理约束：（2.1）考虑阵元间互耦矩阵 (Mutual Coupling)
% 加入后阵元数目为331，副瓣降低至28.7dB,主瓣宽度增至5.0639°，工作空间保存为Mutual_CouplingWorkSpace，迭代轮次多
% （2.2）新增约束：加入干扰的零陷约束
% 3.初始化/网格：
% （3.1）网格：随机扰动网格，开启阵元数目743。阵元数目274，波束宽度5.0639°。工作空间保存为causal_netpoint_WorkSpace，副瓣电平为-27.31dB，迭代29轮
% （3.2）密度锥化阵列初始化：阵元数目261，最高旁瓣电平27.62dB，波束宽度5.0639度，迭代轮次77轮。工作空间保存为Density_coning_WorkSpace.mat
% 4.修改IDEA：改成自适应移动半径：运行结果：阵元数278，最高旁瓣电平27.79dB，波束宽度5.0639度，迭代轮次30轮。工作空间保存为Adaptive_radiusWorkSpace.mat(1.3.4)
% 或者与IDEA相同
% =========================================================================
%1.5范数+随机扰动网格+自适应半径+密度锥化阵列初始化。运行结果：工作空间保存为WorkSpace_all.mat最终阵元数253，迭代35次，副瓣电平-28.42dB，主瓣宽度5.0639度
% =========================================================================
clear; clc; close all;
if exist('cvx_setup.m', 'file'), cvx_clear; end

%% 1. 参数定义
lambda = 1;
beta = 2*pi/lambda;
u1 = 0.12;               % 主瓣半径
SLL_target_dB = -30;    % 【设定】旁瓣目标 -30dB
SLL_limit = 10^(SLL_target_dB/20);
d_min = 0.35 * lambda;   
delta_lock= 0.005;            %锁定点的阵元激励
lambda_L2 = 0.2;        % 能量平滑项系数，防止阵元过快消失
mu_reweight = 1e-3;     % 加权因子的平滑常数
threshold_c = 1e-4;     % 剔除门槛
gamma=0.05;%扰动幅度比
pos_noise_scale = 0.01 * lambda;

% % --- 初始布局：均匀网格 ---
% M_side = 32;            % 初始 22x22 满阵
% d_init = 0.5 * lambda;
% [X_g, Y_g] = meshgrid(-(M_side-1)/2 * d_init : d_init : (M_side-1)/2 * d_init);
% pos_x = X_g(:); pos_y = Y_g(:);
% N = length(pos_x);
% [X_g, Y_g] = meshgrid(1:M_side, 1:M_side); % 记录行列索引
% grid_m = X_g(:); 
% grid_n = Y_g(:);

% --- 初始布局：随机扰动网格 ---
M_side = 32;            
d_init = 0.5 * lambda;
[X_g, Y_g] = meshgrid(-(M_side-1)/2 * d_init : d_init : (M_side-1)/2 * d_init);
% --- 核心修改：注入随机位置抖动 ---
% jitter_scale 决定了偏离规则网格的幅度，通常设为 0.2~0.4 个间距
jitter_scale = 0.3 * d_init; 
pos_x = X_g(:) + (rand(size(X_g(:))) - 0.5) * jitter_scale;
pos_y = Y_g(:) + (rand(size(Y_g(:))) - 0.5) * jitter_scale;
N = length(pos_x);
% 保持逻辑索引记录（即使位置变了，它们在拓扑上依然对应原来的网格点）
[X_idx, Y_idx] = meshgrid(1:M_side, 1:M_side);
grid_m = X_idx(:); 
grid_n = Y_idx(:);



max_iter = 100;
P = 3;                  % 膨胀顶点数
stop_counter = 0;           % 稳定计数器
c =zeros(N, 1) ;   % 初始激励向量

% UV 采样点（高密度采样以保证 -30dB 约束有效）
[U, V] = meshgrid(linspace(-1,1,60));
UV_dist = sqrt(U.^2 + V.^2);
mask_sll = (UV_dist >= u1);
u_s = U(mask_sll); v_s = V(mask_sll);

% %初始化：圆形初始化
% num_elements_history = [];
% R_aperture = (M_side-1)/2* d_init; % 圆形孔径半径
% dist_from_center = sqrt(pos_x.^2 + pos_y.^2);
% in_circle_mask = (dist_from_center <= R_aperture);
% c(in_circle_mask)=1*0.5;
% c=c(in_circle_mask);
% pos_x=pos_x(in_circle_mask);
% pos_y=pos_y(in_circle_mask);
% grid_m = grid_m(in_circle_mask);
% grid_n = grid_n(in_circle_mask);
% fprintf('总网格点数: %d, 开启的圆形阵元数: %d\n', N, sum(in_circle_mask));
% N= sum(in_circle_mask);

%初始化：密度锥化阵列
num_elements_history = [];
R_aperture = (M_side-1)/2 * d_init; % 最大半径
dist_from_center = sqrt(pos_x.^2 + pos_y.^2);
rho = dist_from_center / R_aperture;
% 越靠近中心 (rho->0), 概率越接近 1; 越靠近边缘 (rho->1), 概率越低
n_taper = 0.4; % 锥化指数，值越大，边缘越稀疏
prob_distribution = 0.4 + 0.6 * cos(rho * pi/2).^n_taper;
% prob_distribution(dist_from_center > R_aperture) = 0;% 施加圆形边界掩码：超出半径的概率直接设为 0
density_mask = rand(size(pos_x)) < prob_distribution;
x_max = max(pos_x); x_min = min(pos_x);
y_max = max(pos_y); y_min = min(pos_y);
[~, idx_right] = min(abs(pos_y) + abs(pos_x - x_max));
[~, idx_left]  = min(abs(pos_y) + abs(pos_x - x_min));
[~, idx_top]   = min(abs(pos_x) + abs(pos_y - y_max));
[~, idx_bottom]= min(abs(pos_x) + abs(pos_y - y_min));
lock_idx=[idx_right;idx_left;idx_top;idx_bottom];
density_mask(lock_idx) = true; 
N = sum(density_mask); 
c = 0.5 * ones(N, 1); % 初始激励设为等幅
pos_x = pos_x(density_mask);
pos_y = pos_y(density_mask);
grid_m = grid_m(density_mask);
grid_n = grid_n(density_mask);
fprintf('总网格点数: %d, 密度锥化后初始阵元数: %d (初始稀疏率: %.2f%%)\n', ...
    M_side^2, N, (N/M_side^2)*100);

%% 2. IDEA 核心迭代循环
for k = 1:max_iter
    active_count_old =N ; % 记录上一轮的阵元数
    if N < 10
        fprint('在第%d次迭代开始时阵元数目小于10，停止迭代',k);
        break; 
    end
    %重新定义锁定点
    x_max = max(pos_x); x_min = min(pos_x);
    y_max = max(pos_y); y_min = min(pos_y);
    [~, idx_right] = min(abs(pos_y) + abs(pos_x - x_max));
    [~, idx_left]  = min(abs(pos_y) + abs(pos_x - x_min));
    [~, idx_top]   = min(abs(pos_x) + abs(pos_y - y_max));
    [~, idx_bottom]= min(abs(pos_x) + abs(pos_y - y_min));
    %锁定点的索引
    lock_idx=[idx_right;idx_left;idx_top;idx_bottom];
    %锁定点转化为逻辑值
    lock_mask = false(N,1);
    lock_mask(lock_idx) = true;
    

%     % --- 互耦矩阵动态计算 ---
%     % 或者最小阵元间距约束，使用阵元索引大小排序，严格约束间距大于等于d_min
%     C_mat = zeros(N);
%     coupling_factor=0.1;%耦合强度系数
%     for i = 1:N
%         for j = i+1:N
%             manhattan_dist = abs(grid_m(i) - grid_m(j)) + abs(grid_n(i) - grid_n(j));
%        
%            if manhattan_dist == 1
%             dist_ij = sqrt((pos_x(i)-pos_x(j))^2 + (pos_y(i)-pos_y(j))^2);
%             C_mat(i,j)=dist_ij;
% %             C_mat(i,j) = coupling_factor * exp(-1j*beta*dist_ij)/(dist_ij + 0.1);
% %             C_mat(j,i) = C_mat(i,j);
%             end
%         end
%     end
 if mod(k, 6) == 0
        fprintf('>>> 触发随机扰动：正在为激励和坐标注入抖动...\n');
        % 1. 激励扰动 (对幅度进行 ±5% 左右的随机抖动，并注入微小相位)
        w_noise_amp = gamma * max(abs(c)); 
        c = c + (randn(size(c)) + 1j*randn(size(c))) * (w_noise_amp / 2);
        % 2. 坐标扰动 (注入 ±0.01 lambda 的微小位移）
        pos_x_rand = randn(size(pos_x)) * pos_noise_scale;
        pos_y_rand = randn(size(pos_y)) * pos_noise_scale;

        perturb_mask = true(size(pos_x));
        perturb_mask(lock_idx) = false; % 锁定点不参与坐标抖动
        %加入坐标抖动
        pos_x(perturb_mask) = pos_x(perturb_mask) + pos_x_rand(perturb_mask);
        pos_y(perturb_mask) = pos_y(perturb_mask) + pos_y_rand(perturb_mask);
        %加入抖动后将稳定计数器置零
        stop_counter = 0; 
    end
    
    % --- Step 1: Sparsification (基于 L1.5 范数) ---
    A_sll = exp(1j * beta * (u_s * pos_x' + v_s * pos_y'));
    A_0 = exp(1j * beta * (0 * pos_x' + 0 * pos_y'));
    % --- 稀疏加权向量更新 ---
    % 根据上一次结果 c 更新权重 W。值越小，权重越大，从而促进稀疏
    W_diag = 1 ./ (abs(c) + mu_reweight); 
    
    cvx_clear
    cvx_begin quiet
        variable c_opt(N, 1) complex
        % 目标函数：加权 L1.5 范数 + L2 正则
        minimize( norm(W_diag .* c_opt, 1.5) + lambda_L2 * norm(c_opt, 2) ) 
        subject to
            real(A_0  * c_opt) >= 1; 
            abs(A_sll * c_opt) <= SLL_limit;
            c_opt(lock_idx)==delta_lock;
    cvx_end
    
    % 检查求解状态
    if strcmp(cvx_status, 'Infeasible') || any(isnan(c_opt))
    fprintf('Iter %d: Step 1 不可解，跳过剔除。\n', k);
    end

    keep_idx = abs(c_opt) > threshold_c;
    if sum(keep_idx) < 7
        fprintf('\n[警告] 阵元太少，停止迭代。');
        break; end % 阵元太少则停止
    %更新优化参数
    c = c_opt;

    % --- Step 2: Inflation (膨胀阶段) ---
    N = length(pos_x);
    pos_x_inf = []; pos_y_inf = [];
    map_back = [];   % 记录 inflation → 原始索引
    delta_0 = lambda/60;  
    Gamma_progress=1-(sum(keep_idx)/active_count_old);
    delta_pro=lambda/60;
    delta_min=lambda/60;
    delta_inf=max(delta_min,delta_0+delta_pro*Gamma_progress);
    for n = 1:N
        if lock_mask(n)
            % 锁定点：不膨胀
            pos_x_inf(end+1,1) = pos_x(n);%只有一个元素，（end+1,1)表示列向量
            pos_y_inf(end+1,1) = pos_y(n);
            map_back(end+1,1) = n;%记录位置，也不膨胀
        else
        phi0 = rand * 2*pi;
        for p = 1:P
            idx = (n-1)*P + p;
            pos_x_inf(end+1,1) = pos_x(n) + delta_inf * cos(2*pi*(p-1)/P + phi0);
            pos_y_inf(end+1,1) = pos_y(n) + delta_inf * sin(2*pi*(p-1)/P + phi0);
            map_back(end+1,1) = n;%记录位置，膨胀
        end
        end
    end
    
    % --- Step 3: Inflation Optimization (基于群组 L1.5 范数) ---
    A_inf_sll = exp(1j * beta * (u_s * pos_x_inf' + v_s * pos_y_inf'));
    A_inf_0 = exp(1j * beta * (0 * pos_x_inf' + 0 * pos_y_inf'));
    
    % 膨胀阶段也使用加权。每个阵元的 P 个顶点共享同一个父阵元的权重
    lock_inf_idx = find(lock_mask(map_back));%锁定点索引
    W_group = 1 ./ (abs(c) + mu_reweight);
    cvx_clear
    cvx_begin quiet
        variable c_inf(length(pos_x_inf)) complex
        % 组稀疏目标：每一组的 L1.5 范数之和
        % CVX 中 norms(X, 1.5, 1) 会对每一列求 1.5 范数
        group_l15 = 0;
        for n = 1:N
            idxs = (map_back == n);
            % W_diag(n) 是第 n 个阵元的权重，norm(..., 1.5) 处理该阵元的 P 个膨胀点
            group_l15 = group_l15 + W_group(n) *norm(c_inf(idxs), 1.5);
        end
        minimize( group_l15 + lambda_L2 * norm(c_inf, 2) )
        subject to
            real(A_inf_0 * c_inf) >= 1;
            abs(A_inf_sll * c_inf) <= SLL_limit;
            c_inf(lock_inf_idx) == delta_lock;

    cvx_end

    if strcmp(cvx_status, 'Infeasible')||any(isnan(c_inf))
        fprintf('Iter %d: Step 2 不可解，跳过剔除。\n', k);
        break;
    end

    % --- Step 4: Deflation & Update ---
    new_pos_x = zeros(N, 1); new_pos_y = zeros(N, 1); new_c = zeros(N, 1);
    for n = 1:N
         idxs = find(map_back == n);
        if lock_mask(n)
            % 🔒 锁定点：完全冻结
            new_pos_x(n) = pos_x(n);
            new_pos_y(n) = pos_y(n);
            new_c(n)     = delta_lock;
        else
        v_w = abs(c_inf(idxs));
        if sum(v_w) > threshold_c
            new_pos_x(n) = sum(v_w .* pos_x_inf(idxs)) / sum(v_w);
            new_pos_y(n) = sum(v_w .* pos_y_inf(idxs)) / sum(v_w);
            new_c(n) = sum(c_inf(idxs));
        end
        end
    end
    
    % 剔除逻辑
    active = abs(new_c) > threshold_c;
    active(lock_idx)=true;
    active_count_new = sum(active);
    if sum(active) < 7
        fprintf('\n[警告] 阵元太少，停止迭代。');
        break; % 阵元太少则停止
    end
    pos_x = new_pos_x(active); 
    pos_y = new_pos_y(active); 
    c = new_c(active);
    grid_m=grid_m(active);
    grid_n=grid_n(active);
    N = length(pos_x);
    num_elements_history(k) = sum(active);%动态显示
    fprintf('Iteration %d: Elements = %d, Status = %s\n', k, length(pos_x), cvx_status);
    if active_count_new == active_count_old
        stop_counter = stop_counter + 1;
    else
        stop_counter = 0; % 如果阵元数还在变，重置计数器
    end
    
    % 如果阵元数目连续 4 次保持不变，则认为收敛，跳出循环
    if stop_counter >= 4
        fprintf('迭代收敛，跳出循环: 阵元数稳定在 %d\n', active_count_new);
        break;
    end
end

%% 3. 最终结果评估与绘图
[U_f, V_f] = meshgrid(linspace(-1,1,250));
AF = zeros(size(U_f));
for n = 1:length(pos_x)
    AF = AF + c(n) * exp(1j*beta*(pos_x(n)*U_f + pos_y(n)*V_f));
end
AF_mag = abs(AF)/max(abs(AF(:)));
AF_dB = 20*log10(AF_mag + 1e-6);

figure('Color','w','Position',[100 100 1000 450]);
% 阵列布局图
subplot(1,2,1); 
scatter(pos_x, pos_y, 35, abs(c), 'filled'); 
colormap(jet); colorbar; axis equal; grid on;
xlabel('x/\lambda'); ylabel('y/\lambda');
title(['Sparse Array Layout (N=', num2str(length(pos_x)), ')']);

% 方向图
subplot(1,2,2); 
surf(U_f, V_f, AF_dB, 'EdgeColor', 'none'); 
view(-30, 60); colormap(jet); colorbar;
zlim([-50 0]); xlabel('u'); ylabel('v');
title('Final Pattern (L1.5 Reweighted)');

fprintf('>>> 优化完成。最终阵元数：%d\n', length(pos_x));
%% 计算半功率波束宽度
% --- 提取切片坐标向量 ---
% 这里的 U_f 是 meshgrid 生成的矩阵，我们需要取其一维轴向量
u_sampling = U_f(1, :); 
v_sampling = V_f(:, 1).'; % 转置为行向量以保持一致

% --- 提取中心切片 (假设主瓣在 u=0, v=0) 
[~, center_idx_u] = min(abs(u_sampling - 0));
[~, center_idx_v] = min(abs(v_sampling - 0));

% 提取并归一化切片数据，防止因最大值不是0dB导致函数无法找到-3dB点
slice_u_dB = AF_dB(center_idx_v, :);
slice_u_dB = slice_u_dB - max(slice_u_dB); % 强制归一化

slice_v_dB = AF_dB(:, center_idx_u).'; % 提取列并转置为向量
slice_v_dB = slice_v_dB - max(slice_v_dB); % 强制归一化

% --- 调用函数计算 HPBW ---
hpbw_u = calculate_hpbw(u_sampling, slice_u_dB);
hpbw_v = calculate_hpbw(v_sampling, slice_v_dB);

% u轴旁瓣标注逻辑 (基于 dB 数据寻峰)
[all_pks, all_locs] = findpeaks(slice_u_dB, u_sampling);
% 识别并剔除主瓣 (通常是最大的峰值)
[~, mainpeak_idx] = max(all_pks); 
if ~isempty(mainpeak_idx)
    all_pks(mainpeak_idx) = -inf; 
end
% 寻找最高旁瓣
[max_sll_val, max_sll_idx] = max(all_pks);
u_max_sll = all_locs(max_sll_idx);

% v轴旁瓣标注逻辑 (基于 dB 数据寻峰)
[all_pks_v, all_locs_v] = findpeaks(slice_v_dB, v_sampling);
% 识别并剔除主瓣 (峰值最大处)
[~, mainpeak_idx_v] = max(all_pks_v); 
if ~isempty(mainpeak_idx_v)
    all_pks_v(mainpeak_idx_v) = -inf; 
end
% 寻找最高旁瓣
[max_v_sll_val, max_v_sll_idx] = max(all_pks_v);
v_max_sll_pos = all_locs_v(max_v_sll_idx);


% --- 格式化输出结果 ---
fprintf('\n================ 阵列性能评估指标 =================\n');
fprintf('最终有效阵元数目 (N):   %d\n', length(pos_x));
fprintf('稀疏率 (Sparse Ratio):  %.2f%%\n', (length(pos_x)/M_side^2)*100);
fprintf('3dB 波束宽度 (u-axis):  %.4f 度\n', hpbw_u);
fprintf('3dB 波束宽度 (v-axis):  %.4f 度\n', hpbw_v);
fprintf('最高副瓣电平 (u-axis):  %.2f dB\n', max_sll_val);
fprintf('最高副瓣电平 (v-axis):  %.2f dB\n', max_v_sll_val);
fprintf('===================================================\n');



%% 相同副瓣的满阵加窗比较
 %u方向函数引用
[U_f_taylor,U_f_ChebyShev] = full_array_with_windows_aftersamesidelobe_u(M_side,beta,d_init,max_sll_val,250);
 %v方向函数引用
[V_f_taylor,V_f_ChebyShev] = full_array_with_windows_aftersamesidelobe_v(M_side,beta,d_init,max_v_sll_val,250);

%u轴切片方向图显示
figure;
plot(u_sampling, slice_u_dB, 'LineWidth', 1.8,'DisplayName', '稀疏优化阵列');
hold on; % 
plot(u_sampling, U_f_taylor, 'm--', 'LineWidth', 1.8, 'DisplayName', '满阵加泰勒窗');
plot(u_sampling, U_f_ChebyShev, 'c--', 'LineWidth', 1.8, 'DisplayName', '满阵加切比雪夫窗');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴基本属性
grid on;
ylim([-40 5]); % 建议上限给到 5dB，防止主瓣顶端被压住
xlim([min(u_sampling) max(u_sampling)]);
xlabel('u = sin(\theta)cos(\phi)');
ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）（相同旁瓣电平）');
hold off;

%v轴切片方向图显示
figure;
% 使用橙红色系曲线以区分 u 轴切片
plot(v_sampling, slice_v_dB, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980],'DisplayName', '稀疏优化阵列'); 
hold on;
plot(v_sampling, V_f_taylor, 'b--', 'LineWidth', 1.8, 'DisplayName', '满阵加泰勒窗');
plot(v_sampling, V_f_ChebyShev, 'y--', 'LineWidth', 1.8, 'DisplayName', '满阵加切比雪夫窗');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴属性
grid on;
ylim([-40 5]); 
xlim([min(v_sampling) max(v_sampling)]);
xlabel('v = sin(\theta)sin(\phi)');
ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)（相同旁瓣电平）');

%分别计算3dB波束宽度
[hpbw_sparse_u] = calculate_hpbw(u_sampling, slice_u_dB);
[hpbw_taylor_u] = calculate_hpbw(u_sampling, U_f_taylor);
[hpbw_ChebyShev_u] = calculate_hpbw(u_sampling,U_f_ChebyShev);

fprintf("----------------相同副瓣的满阵加窗比较------------------------\n");
fprintf('u轴稀疏阵列 (3dB宽度): %.4f°\n', hpbw_sparse_u);
fprintf('u轴满阵加泰勒窗 (3dB宽度): %.4f°\n', hpbw_taylor_u);
fprintf('u轴满阵加切比雪夫窗 (3dB宽度): %.4f°\n', hpbw_ChebyShev_u);

[hpbw_sparse_v] = calculate_hpbw(v_sampling, slice_v_dB);
[hpbw_taylor_v] = calculate_hpbw(v_sampling, V_f_taylor);
[hpbw_ChebyShev_v] = calculate_hpbw(v_sampling,V_f_ChebyShev);

fprintf('v轴稀疏阵列 (3dB宽度): %.4f°\n', hpbw_sparse_v);
fprintf('v轴满阵加泰勒窗 (3dB宽度): %.4f°\n', hpbw_taylor_v);
fprintf('v轴满阵加切比雪夫窗 (3dB宽度): %.4f°\n', hpbw_ChebyShev_v);
%% 其他方法的对比
%--------------------------------------------------------------------------
[U_f, V_f] = meshgrid(linspace(-1,1,1000));
AF = zeros(size(U_f));
for n = 1:length(pos_x)
    AF = AF + c(n) * exp(1j*beta*(pos_x(n)*U_f + pos_y(n)*V_f));
end
AF_mag = abs(AF)/max(abs(AF(:)));
AF_dB = 20*log10(AF_mag + 1e-6);
u_sampling = U_f(1, :); 
v_sampling = V_f(:, 1).'; % 转置为行向量以保持一致

% --- 提取中心切片 (假设主瓣在 u=0, v=0) 
[~, center_idx_u] = min(abs(u_sampling - 0));
[~, center_idx_v] = min(abs(v_sampling - 0));

% 提取并归一化切片数据，防止因最大值不是0dB导致函数无法找到-3dB点
slice_u_dB = AF_dB(center_idx_v, :);
slice_u_dB = slice_u_dB - max(slice_u_dB); % 强制归一化

slice_v_dB = AF_dB(:, center_idx_u).'; % 提取列并转置为向量
slice_v_dB = slice_v_dB - max(slice_v_dB); % 强制归一化
%--------------------------------------------------------------------------
%IDEA
folderpath='C:\Users\gyq\Desktop\sparse_arrays_design\others_method\IDEA';
filename='IDEA5WorkSpace.mat';
fullPath=fullfile(folderpath, filename);

dataStruct=load(fullPath,'U_f','V_f','pos_x','pos_y','AF_dB','AF','c','slice_axis','pattern_u','pattern_v','show_sampling');

fprintf('IDEA方法最终阵元数目为：%d\n',length(dataStruct.pos_x));

%u轴切片方向图显示
figure;
plot(u_sampling, slice_u_dB, 'LineWidth', 1.8,'DisplayName', '稀疏优化阵列');
hold on; % ！！！关键点：必须 hold on 才能在同一张图上画标注
plot(u_sampling,dataStruct.pattern_u, 'm--', 'LineWidth', 1.8, 'DisplayName', 'IDEA');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on;
ylim([-40 5]); % 建议上限给到 5dB，防止主瓣顶端被压住
xlim([min(u_sampling) max(u_sampling)]);
xlabel('u = sin(\theta)cos(\phi)');
ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）（不同优化方法比较）');
hold off;
%v轴切片方向图显示
figure;
plot(v_sampling, slice_v_dB, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980],'DisplayName', '稀疏优化阵列'); % 使用橙红色系曲线以区分 u 轴切片
hold on;
plot(v_sampling, dataStruct.pattern_v, 'b--', 'LineWidth', 1.8, 'DisplayName', 'IDEA');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on;
ylim([-40 5]); 
xlim([min(v_sampling) max(v_sampling)]);
xlabel('v = sin(\theta)sin(\phi)');
ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)（不同优化方法比较）');

%% %% 和其他优化方法比较主瓣增益
% U_power = abs(AF).^2;
% U_max = max(U_power(:));
% du = 2/(1000-1);
% dv = 2/(1000-1);
% valid_mask = (U_f.^2 + V_f.^2) < 1; % 建议用 < 1 而非 <=1，防止分母为0导致 Inf
% cos_theta = sqrt(1 - U_f(valid_mask).^2 - V_f(valid_mask).^2);
% P_rad = sum( U_power(valid_mask) ./ cos_theta ) * (du * dv); 
% % 计算最终的方向性 (此时才等效于严格的球面积分)
% Directivity_linear = (4 * pi * U_max) / P_rad;
% Directivity_dB = 10 * log10(Directivity_linear);
% fprintf('本稀疏阵列优化方法的严谨方向性系数为: %.2f dB\n', Directivity_dB);
% 
% 
% U_power_IDEA = abs(dataStruct.AF).^2;
% U_max_IDEA = max(U_power_IDEA(:));
% du_IDEA = 2/(dataStruct.show_sampling-1);
% dv_IDEA = 2/(dataStruct.show_sampling-1);
% valid_mask_IDEA = (dataStruct.U_f.^2 + dataStruct.V_f.^2) <= 1;
% cos_theta_IDEA = sqrt(1 - dataStruct.U_f(valid_mask_IDEA).^2 - dataStruct.V_f(valid_mask_IDEA).^2);
% P_rad_IDEA = sum(U_power_IDEA(valid_mask_IDEA)./cos_theta_IDEA) * (du_IDEA * dv_IDEA); % du, dv 为采样步长
% Directivity_linear_IDEA = (4 * pi * U_max_IDEA) / P_rad_IDEA;
% Directivity_dB_IDEA = 10 * log10(Directivity_linear_IDEA);
% fprintf('IDEA优化方法的方向性系数为: %.2f dB\n', Directivity_dB_IDEA);
%% %% 和其他优化方法比较主瓣增益 (切换至严谨的球坐标系积分)

% --- 参数设置 ---
d_theta = 0.5; % theta 采样步长（度）
d_phi = 1.0;   % phi 采样步长（度）
% 积分范围：theta从0到90度(上半球面)，phi从0到359度
[TH, PH] = meshgrid(deg2rad(0:d_theta:90), deg2rad(0:d_phi:359));

% 映射为 u, v 用于计算阵列因子
U_sph = sin(TH) .* cos(PH);
V_sph = sin(TH) .* sin(PH);

% ---------------------------------------------------------
% 1. 计算本稀疏阵列的方向性系数
% ---------------------------------------------------------
% 强制转换为列向量以使用矩阵相乘 (提速)
pos_x_col = pos_x(:);
pos_y_col = pos_y(:);
c_col = c(:);

% 矩阵相乘快速计算球面所有点的 AF: (1xN) * (NxM)
AF_sph = abs(c_col.' * exp(1j * beta * (pos_x_col * U_sph(:).' + pos_y_col * V_sph(:).')));
U_power_sph = AF_sph.^2;

U_max = max(U_power_sph(:));

% 球面积分: sum(|AF|^2 * sin(theta) * d_theta * d_phi)
sin_theta_vec = sin(TH(:)).'; 
P_rad = sum(U_power_sph .* sin_theta_vec) * deg2rad(d_theta) * deg2rad(d_phi);

Directivity_linear = (4 * pi * U_max) / P_rad;
Directivity_dB = 10 * log10(Directivity_linear);
fprintf('本稀疏阵列优化方法 (球坐标积分) 的方向性系数为: %.2f dBi\n', Directivity_dB);


% ---------------------------------------------------------
% 2. 计算 IDEA 优化方法的方向性系数
% ---------------------------------------------------------
if exist('dataStruct', 'var')
    pos_x_IDEA = dataStruct.pos_x(:);
    pos_y_IDEA = dataStruct.pos_y(:);
    c_IDEA = dataStruct.c(:);

    AF_sph_IDEA = abs(c_IDEA.' * exp(1j * beta * (pos_x_IDEA * U_sph(:).' + pos_y_IDEA * V_sph(:).')));
    U_power_sph_IDEA = AF_sph_IDEA.^2;

    U_max_IDEA = max(U_power_sph_IDEA(:));

    P_rad_IDEA = sum(U_power_sph_IDEA .* sin_theta_vec) * deg2rad(d_theta) * deg2rad(d_phi);

    Directivity_linear_IDEA = (4 * pi * U_max_IDEA) / P_rad_IDEA;
    Directivity_dB_IDEA = 10 * log10(Directivity_linear_IDEA);
    fprintf('IDEA优化方法 (球坐标积分) 的方向性系数为: %.2f dBi\n', Directivity_dB_IDEA);
end
%% 和满阵进行比较
%--------------------------------------------------------------------------
% IDEA 方法数据载入
folderpath='C:\Users\gyq\Desktop\sparse_arrays_design\others_method\IDEA';
filename='IDEA5WorkSpace.mat';
fullPath=fullfile(folderpath, filename);
dataStruct=load(fullPath,'U_f','V_f','pos_x','pos_y','AF_dB','AF','c','slice_axis','pattern_u','pattern_v','show_sampling');
fprintf('IDEA方法最终阵元数目为：%d\n',length(dataStruct.pos_x));

%--------------------------------------------------------------------------
% 计算 32x32 均匀满阵的切片数据
pos_x_full = X_g(:);
pos_y_full = Y_g(:);
% u方向 (v=0)
AF_full_u_complex = sum(exp(1j*beta*(pos_x_full * u_sampling)), 1);
AF_full_u_dB = 20*log10(abs(AF_full_u_complex) / max(abs(AF_full_u_complex)));
% v方向 (u=0)
AF_full_v_complex = sum(exp(1j*beta*(pos_y_full * v_sampling)), 1);
AF_full_v_dB = 20*log10(abs(AF_full_v_complex) / max(abs(AF_full_v_complex)));

%--------------------------------------------------------------------------
% 统一定义三种主题色
color_sparse = [0.8500, 0.3250, 0.0980]; % 橙红色 (本方法)
color_idea   = [0.4660, 0.6740, 0.1880]; % 绿色 (IDEA方法)
color_full   = [0.0000, 0.4470, 0.7410]; % 经典蓝 (满阵)

% ========== 1. U 轴切片方向图对比与自动标注 ==========
figure('Color', 'w', 'Position', [100 100 800 500]);
% 绘制三条曲线
plot(u_sampling, slice_u_dB, 'LineWidth', 1.8, 'Color', color_sparse, 'DisplayName', '本稀疏优化阵列');
hold on;
plot(u_sampling, dataStruct.pattern_u, '-.', 'LineWidth', 1.8, 'Color', color_idea, 'DisplayName', 'IDEA 方法');
plot(u_sampling, AF_full_u_dB, '--', 'LineWidth', 1.5, 'Color', color_full, 'DisplayName', '均匀满阵');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on; ylim([-45 5]); xlim([min(u_sampling) max(u_sampling)]);
xlabel('u = sin(\theta)cos(\phi)'); ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）（多方法对比）');

% --- 自动寻峰与标注模块 (U轴) ---
data_list_u = {slice_u_dB, dataStruct.pattern_u, AF_full_u_dB};
color_list = {color_sparse, color_idea, color_full};
name_list = {'稀疏阵', 'IDEA', '满阵'};
y_offset_list = [1.5, -2.5, 1.5]; % 分别为三种文字设置 Y 轴偏移，防止互相重叠

for i = 1:3
    [pks, locs] = findpeaks(data_list_u{i}, u_sampling);
    [~, main_idx] = max(pks); 
    if ~isempty(main_idx), pks(main_idx) = -inf; end % 剔除主瓣
    [max_sll, sll_idx] = max(pks);
    sll_pos = locs(sll_idx);
    
    if ~isempty(max_sll)
        % 画圆圈标记
        plot(sll_pos, max_sll, 'o', 'MarkerEdgeColor', color_list{i}, 'MarkerSize', 8, 'LineWidth', 2, 'HandleVisibility', 'off');
        % 画参考线
        line([min(u_sampling) max(u_sampling)], [max_sll max_sll], 'Color', color_list{i}, 'LineStyle', ':', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        % 动态调整文字左右对齐
        if sll_pos > 0, align = 'left'; x_off = 0.03; else, align = 'right'; x_off = -0.03; end
        % 打印文字
        text(sll_pos + x_off, max_sll + y_offset_list(i), ...
            sprintf('%s: %.2f dB', name_list{i}, max_sll), ...
            'Color', color_list{i}, 'FontWeight', 'bold', 'HorizontalAlignment', align);
    end
end
hold off;

% ========== 2. V 轴切片方向图对比与自动标注 ==========
figure('Color', 'w', 'Position', [150 150 800 500]);
% 绘制三条曲线
plot(v_sampling, slice_v_dB, 'LineWidth', 1.8, 'Color', color_sparse, 'DisplayName', '本稀疏优化阵列'); 
hold on;
plot(v_sampling, dataStruct.pattern_v, '-.', 'LineWidth', 1.8, 'Color', color_idea, 'DisplayName', 'IDEA 方法');
plot(v_sampling, AF_full_v_dB, '--', 'LineWidth', 1.5, 'Color', color_full, 'DisplayName', '均匀满阵');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on; ylim([-45 5]); xlim([min(v_sampling) max(v_sampling)]);
xlabel('v = sin(\theta)sin(\phi)'); ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)（多方法对比）');

% --- 自动寻峰与标注模块 (V轴) ---
data_list_v = {slice_v_dB, dataStruct.pattern_v, AF_full_v_dB};

for i = 1:3
    [pks, locs] = findpeaks(data_list_v{i}, v_sampling);
    [~, main_idx] = max(pks); 
    if ~isempty(main_idx), pks(main_idx) = -inf; end % 剔除主瓣
    [max_sll, sll_idx] = max(pks);
    sll_pos = locs(sll_idx);
    
    if ~isempty(max_sll)
        % 画圆圈标记
        plot(sll_pos, max_sll, 'o', 'MarkerEdgeColor', color_list{i}, 'MarkerSize', 8, 'LineWidth', 2, 'HandleVisibility', 'off');
        % 画参考线
        line([min(v_sampling) max(v_sampling)], [max_sll max_sll], 'Color', color_list{i}, 'LineStyle', ':', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        % 动态调整文字左右对齐
        if sll_pos > 0, align = 'left'; x_off = 0.03; else, align = 'right'; x_off = -0.03; end
        % 打印文字
        text(sll_pos + x_off, max_sll + y_offset_list(i), ...
            sprintf('%s: %.2f dB', name_list{i}, max_sll), ...
            'Color', color_list{i}, 'FontWeight', 'bold', 'HorizontalAlignment', align);
    end
end
hold off;

