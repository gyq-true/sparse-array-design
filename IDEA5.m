%% IDEA Algorithm for Sparse Planar Array Synthesis - Robust Version
%大规模版本
%问题：不严格约束、未加入最小阵元间距约束、忘记加锁定点约束了
%可能的原因：采样间隔过大
%措施：
clear; clc; close all;
if exist('cvx_setup.m', 'file'), cvx_clear; end

%% 1. 更加宽松的初始参数
lambda = 1;
beta = 2*pi/lambda;
u1 = 0.1;               % 【增大主瓣宽度】初始设为 0.2，更容易找到可行解
SLL_target_dB = -30;    % 【放宽旁瓣要求】先从 -15dB 开始尝试，跑通后再调优
SLL_limit = 10^(SLL_target_dB/20);
delta_lock= 0.005;            %锁定点的阵元激励
P = 3; % 膨胀多边形顶点数 (平面阵列选三角形)

% --- 新增：最小间距约束处理 ---
d_min = 0.5 * lambda; % 设定最小间距阈值


% 初始网格：稍微加大一点初始规模，确保有足够的自由度
M_side = 32;            % 改为 10x10 = 100 阵元
d_init = 0.5 * lambda;
[X, Y] = meshgrid(-(M_side-1)/2 * d_init : d_init : (M_side-1)/2 * d_init, ...
                  -(M_side-1)/2 * d_init : d_init : (M_side-1)/2 * d_init);
pos_x = X(:); pos_y = Y(:);
N = length(pos_x);

% 优化控制
max_iter = 10;
delta = lambda/60;
mu = 1e-4;
epsilon = 0.05;        % 增大剔除阈值，加快稀疏速度
threshold=1e-8;%权重剔除阈值

% 旁瓣区域采样 (减小采样点数提高 CVX 成功率)
[U, V] = meshgrid(linspace(-1,1,50), linspace(-1,1,50));
UV_dist = sqrt(U.^2 + V.^2);
mask=(UV_dist >= u1);
%mask = (UV_dist >= u1) & (UV_dist <= 1);%环形旁瓣区域
u_s = U(mask); v_s = V(mask);%旁瓣区域


% % --- 可视化 UV 采样点分布 ---
% figure('Color', 'w', 'Name', 'UV Sampling Points Distribution');
% hold on;
% 
% % 1. 绘制单位圆边界 (初始化圆形阵列)
% theta_circle = linspace(0, 2*pi, 100);
% plot(cos(theta_circle), sin(theta_circle), 'k--', 'LineWidth', 1.5);
% 
% % 2. 绘制主瓣约束区域 (半径为 u1)
% plot(u1*cos(theta_circle), u1*sin(theta_circle), 'r', 'LineWidth', 2);
% 
% % 3. 绘制采样点 (u_s, v_s)
% scatter(u_s, v_s, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
% 
% % 4. 设置坐标轴和标注
% xlabel('u = sin\theta cos\phi');
% ylabel('v = sin\theta sin\phi');
% title(['UV 可视化分布 (Num = ', num2str(length(u_s)), ')']);
% legend('初始化圆形阵列 (u^2+v^2=1)', 'Mainlobe Region (u < u1)', 'Sidelobe Sampling Points');
% axis equal;
% grid on;
% xlim([-1.1 1.1]); ylim([-1.1 1.1]);
% hold off;

num_elements_history = [];
c = zeros(N, 1); %初始阵元激励

%通过切比雪夫阵列计算起始阵元数目
% M=1+(acosh(1/abs(SLL_target_dB))/acosh(1/cos(u1*pi/2)));
% M_R=sqrt(M/4*pi);
% --- 修改点：网格生成后的圆形裁剪 ---
%初始化圆形区域内阵元激励均为1，锁定点阵元激励也为1（包括了）
R_aperture = (M_side-1)/2* d_init; % 圆形孔径半径
dist_from_center = sqrt(pos_x.^2 + pos_y.^2);
in_circle_mask = (dist_from_center <= R_aperture);
c(in_circle_mask)=1;
fprintf('总网格点数: %d, 开启的圆形阵元数: %d\n', N, sum(in_circle_mask));

%找到主瓣中心的点
% 找到距离最小的4个点的索引 (对于偶数网格，这4个点距离相等)
[~, sorted_idx] = sort(dist_from_center);
origin_neighbors_idx = sorted_idx(1:4); 

%% 2. 迭代优化
for k = 1:max_iter
    cvx_clear;
    %找到四个孔径保持点的索引，这里选择区域和坐标轴的四个交点
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
    free_mask = ~lock_mask;
    
    %循环内重新定义主瓣中心点
    dist_of_sparse = sqrt(pos_x.^2 + pos_y.^2);
    [~, sorted_idx] = sort(dist_of_sparse);
    origin_neighbors_idx = sorted_idx(1:4); 
    %% 第一个凸优化
    N=length(pos_x);
    A_sll = exp(1j * beta * (u_s * pos_x' + v_s * pos_y'));
    A_0 = exp(1j * beta * (0* pos_x'+ 0 * pos_y'));
    W = 1 ./ (abs(c) + mu);    %加权
    %检查锁定点与主瓣约束是否冲突
    if any(ismember(origin_neighbors_idx, lock_idx))
    error('主瓣约束点与锁定点发生冲突');
    end

    cvx_begin quiet
        variable c_opt(N,1) complex
        minimize( norm(W.* c_opt, 1) )
        subject to
            real(A_0*c_opt)>= 1; 
            abs(A_sll * c_opt) <= SLL_limit;
            c_opt(lock_idx)==delta_lock;
    cvx_end
    
    % --- 鲁棒性检查 ---
    if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
        fprintf('\n[警告] 迭代 %d 求解失败。', k);
        %fprintf('\n建议：增大 u1 或提高 SLL_target_dB。\n');
        break; % 停止迭代，显示已有结果
    end
    
    % 更新阵元激励
   % keep_idx = abs(c_opt) > epsilon * max(abs(c_opt));
    keep_idx = abs(c_opt) > threshold;
    if sum(keep_idx) < 7
        fprintf('\n[警告] 阵元太少，停止迭代。');
        break; end % 阵元太少则停止
    %更新优化参数
    c = c_opt;

    %最小阵元间距约束
%     N_temp = length(pos_x);
%     to_remove = false(N_temp, 1);
%     for i = 1:N_temp
%         if to_remove(i), continue; end
%         for j = i+1:N_temp
%             % 计算第 i 和第 j 个阵元之间的距离
%             dist = sqrt((pos_x(i)-pos_x(j))^2 + (pos_y(i)-pos_y(j))^2);
%             
%             if dist < d_min
%                 % 如果靠得太近，标记励磁较小的那个为“待剔除”
%                 if abs(c(i)) >= abs(c(j))
%                     to_remove(j) = true;
%                 else
%                     to_remove(i) = true;
%                     break; % i 已被剔除，跳出内层循环
%                 end
%             end
%         end
%     end
    
    % 执行剔除
%     pos_x(to_remove) = [];
%     pos_y(to_remove) = [];
%     c(to_remove) = []; 
%     % -----------------------
  
    %% --- Step 2: Inflation (阵元膨胀与探测) ---
    % 此时 pos_x, pos_y 包含上一轮的 N个阵元位置
    N = length(pos_x);
   
    %阵元膨胀（除去锁定点）
    pos_x_inf = [];
    pos_y_inf = [];
    map_back = [];   % 记录 inflation → 原始索引
    
    for n = 1:N
        if lock_mask(n)
            % 锁定点：不膨胀
            pos_x_inf(end+1,1) = pos_x(n);%只有一个元素，（end+1,1)表示列向量
            pos_y_inf(end+1,1) = pos_y(n);
            map_back(end+1,1) = n;%记录位置，也不膨胀
        else
            % 自由点：正常膨胀
            phi0 = rand * 2*pi; % 随机初始旋转角，增加搜索多样性
            for p = 1:P
                pos_x_inf(end+1,1) = pos_x(n) + delta*cos(2*pi*(p-1)/P+phi0);
                pos_y_inf(end+1,1) = pos_y(n) + delta*sin(2*pi*(p-1)/P+phi0);
                map_back(end+1,1) = n;%记录位置，膨胀

            end
        end
    end
    
    % 1. 生成膨胀后的顶点坐标 (维度变为 N * P)
%     pos_x_inf = zeros(N * P, 1);
%     pos_y_inf = zeros(N * P, 1);
%     
    
%     for n = 1:N
%         phi0 = rand * 2*pi; % 随机初始旋转角，增加搜索多样性
%         for p = 1:P
%             idx = (n-1)*P + p;
%             theta_p = 2*pi*(p-1)/P + phi0;
%             pos_x_inf(idx) = pos_x(n) + delta * cos(theta_p);
%             pos_y_inf(idx) = pos_y(n) + delta * sin(theta_p);
%         end
%     end

    lock_inf_idx = find(lock_mask(map_back));%锁定点索引
    center_inf_idx = find(ismember(map_back, origin_neighbors_idx));%主瓣等效点索引

    % 2. 在膨胀后的坐标上重新进行凸优化 (寻找最优探测点)
    A_inf = exp(1j * beta * (u_s * pos_x_inf' + v_s * pos_y_inf'));
    A_inf_0=exp(1j * beta * (0 * pos_x_inf' +  0 * pos_y_inf'));
    % 使用 kron 函数将 c 的每个元素扩展为 P 个相同的元素并平分幅值
    %c_inf_init = kron(c, ones(P, 1) )/ P;
    % 修改 c_inf_init 的生成逻辑
    c_inf_init = zeros(length(pos_x_inf), 1);
    for n = 1:N
        % 找到膨胀空间中属于原始第 n 个点的所有索引
        current_inf_indices = (map_back == n); 
        num_sub_elements = sum(current_inf_indices);
        % 平分原始激励
        c_inf_init(current_inf_indices) = c(n) / num_sub_elements;
    end
    
    W_inf = 1 ./ (abs(c_inf_init) + mu);

    cvx_begin quiet
    variable c_inf(length(pos_x_inf)) complex
    minimize(norm(W_inf.* c_inf, 1))
    subject to
        real(A_inf_0*c_inf)>= 1;
        abs(A_inf * c_inf) <= SLL_limit;
        c_inf(lock_inf_idx) == delta_lock;
    cvx_end
        
    % 3. Deflation (收缩/位置更新): 核心重心法
    new_pos_x = pos_x;
    new_pos_y = pos_y;
    new_c     = zeros(N,1);
    
    for n = 1:N
        idxs = find(map_back == n);
        if lock_mask(n)
            % 🔒 锁定点：完全冻结
            new_pos_x(n) = pos_x(n);
            new_pos_y(n) = pos_y(n);
            new_c(n)     = delta_lock;
        else
            weights = abs(c_inf(idxs));
            if sum(weights) > threshold% 避免除以 0, (如果保持量级的话应该c_inf*3）（保持一致的话应该使用动态的方法）
                % 加权形心更新位置 (Deflation)
                new_pos_x(n) = sum(weights .* pos_x_inf(idxs)) / sum(weights);
                new_pos_y(n) = sum(weights .* pos_y_inf(idxs)) / sum(weights);
                new_c(n)     = sum(c_inf(idxs));
            else
            % 如果该组顶点全为 0，标记为待剔除
            new_pos_x(n) = pos_x(n);
            new_pos_y(n) = pos_y(n);
            new_c(n) = 0;
            end
        end
    end
    
    % 4. 更新下一轮迭代的输入（此处按照激励更新权重的长度）（迭代末尾更新长度）
    % 定义一个活跃掩码：只有激励幅值大于阈值的阵元才保留
    %active_mask = (abs(new_c) > epsilon * max(abs(new_c))) | lock_mask;%动态剔除
    active_mask = abs(new_c) > threshold ;
    active_mask(lock_idx)=true;
    if sum(active_mask) < 7
        fprintf('\n[警告] 阵元太少，停止迭代。');
        break; % 阵元太少则停止
    end
    
    % 更新位置和激励，数组长度在这里会发生收缩
    pos_x = new_pos_x(active_mask);
    pos_y = new_pos_y(active_mask);
    c = new_c(active_mask);% 更新后的 c 将作为下一轮 Sparsification 的权重依据 W
        
    %阵元激励更新
    N = length(pos_x);
    %记录每个迭代的活跃阵元数目
    num_elements_history(k) = sum(active_mask);%动态显示
    
    fprintf('Iteration %d: Elements = %d, Status = %s\n', k, num_elements_history(k), cvx_status);
end

%% 3. 增强版可视化逻辑
% --- 准备数据 ---
show_sampling=200;
[U_f, V_f] = meshgrid(linspace(-1,1,show_sampling), linspace(-1,1,show_sampling)); % 提高分辨率
AF = zeros(size(U_f));
for n = 1:length(pos_x)
    AF = AF + c(n) * exp(1j*beta*(pos_x(n)*U_f + pos_y(n)*V_f));
end
AF_dB = 20*log10(abs(AF)/max(abs(AF(:))));
AF_dB(AF_dB < -60) = -60; % 限制底噪，增强视觉效果

% =========================================================================
% FIGURE 1: 阵元分布与三维方向图
% =========================================================================
figure('Color','w','Name','Array Layout and 3D Pattern','Position',[100 100 1100 500]);

% 左侧：阵元分布
% subplot(1,2,1);
% scatter(pos_x, pos_y, 40, abs(c), 'filled');
% % === 用红叉标出锁定阵元 ===
% scatter(pos_x(lock_mask), pos_y(lock_mask), ...
%         120, 'r', 'x', 'LineWidth', 2);
% xlabel('x / \lambda'); ylabel('y / \lambda');
% title(['Final Sparse Layout (N = ', num2str(length(pos_x)), ')']);
% colorbar; axis equal; grid on;
% 左侧：阵元分布（相对激励强度）
subplot(1,2,1);
c_abs = abs(c);                 % 阵元激励绝对值
c_min = min(c_abs);
c_max = max(c_abs);
% 归一化到 [0,1]
c_rel = (c_abs - c_min) / (c_max - c_min + eps);
scatter(pos_x, pos_y, 40, c_rel, 'filled');
% === 用红叉标出锁定阵元 ===
hold on;
x_max = max(pos_x); x_min = min(pos_x);
y_max = max(pos_y); y_min = min(pos_y);
[~, idx_right] = min(abs(pos_y) + abs(pos_x - x_max));
[~, idx_left]  = min(abs(pos_y) + abs(pos_x - x_min));
[~, idx_top]   = min(abs(pos_x) + abs(pos_y - y_max));
[~, idx_bottom]= min(abs(pos_x) + abs(pos_y - y_min));
%锁定点的索引
lock_idx=[idx_right;idx_left;idx_top;idx_bottom];
scatter(pos_x(lock_idx), pos_y(lock_idx), ...
        120, 'r', 'x', 'LineWidth', 2);
hold off;
xlabel('x / \lambda');
ylabel('y / \lambda');
title(['Final Sparse Layout (N = ', num2str(length(pos_x)), ')']);
colormap(jet);
cb = colorbar;
cb.Label.String = 'Relative Excitation |c| (min \rightarrow max)';
axis equal;
grid on;


% 右侧：3D 方向图
subplot(1,2,2);
surf(U_f, V_f, AF_dB, 'EdgeColor', 'none');
xlabel('u = sin\theta cos\phi'); 
ylabel('v = sin\theta sin\phi');
zlabel('Pattern Magnitude (dB)');
title('3D Radiation Pattern');
view(-35, 45); colorbar; colormap jet;

% =========================================================================
% FIGURE 2: 迭代过程分析
% =========================================================================
figure('Color','w','Name','Iteration Process','Position',[200 200 600 450]);
plot(1:length(num_elements_history), num_elements_history, '-ks', ...
    'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('Iteration Number (k)'); 
ylabel('Number of Active Elements (N_{active})');
title('Convergence: Element Count vs. Iterations');
grid on; set(gca, 'FontSize', 10);

% =========================================================================
% FIGURE 3: U轴与V轴切片分析
% =========================================================================
figure('Color','w','Name','Pattern Cut Analysis','Position',[300 300 1000 500]);
slice_axis = linspace(-1,1,1000);

% (1) U-axis cut (at v=0)
subplot(1,2,1);
A_u_slice = exp(1j * beta * (slice_axis' * pos_x'));
pattern_u = 20*log10(abs(A_u_slice * c)/max(abs(A_u_slice * c)));
plot(slice_axis, pattern_u, 'b', 'LineWidth', 1.5); hold on;
yline(SLL_target_dB, '--r', 'Threshold', 'LineWidth', 1.2);
xlabel('u (sin\theta when \phi=0)'); ylabel('Magnitude (dB)');
title('U-Axis Cut (v=0)'); grid on; ylim([-60 5]);

% 标注 U 轴 PSL
[pks, locs] = findpeaks(pattern_u, slice_axis);
s_mask = abs(locs) > u1;
if any(s_mask)
    [psl, idx] = max(pks(s_mask));
    s_locs = locs(s_mask);
    plot(s_locs(idx), psl, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    text(s_locs(idx)+0.05, psl+2, sprintf('PSL_{u}=%.1f dB', psl), 'Color', 'r', 'FontWeight', 'bold');
end

% (2) V-axis cut (at u=0)
subplot(1,2,2);
A_v_slice = exp(1j * beta * (slice_axis' * pos_y'));
pattern_v = 20*log10(abs(A_v_slice * c)/max(abs(A_v_slice * c)));
plot(slice_axis, pattern_v, 'm', 'LineWidth', 1.5); hold on;
yline(SLL_target_dB, '--r', 'Threshold', 'LineWidth', 1.2);
xlabel('v (sin\theta when \phi=90^\circ)'); ylabel('Magnitude (dB)');
title('V-Axis Cut (u=0)'); grid on; ylim([-60 5]);

% 标注 V 轴 PSL
[pks, locs] = findpeaks(pattern_v, slice_axis);
s_mask = abs(locs) > u1;
if any(s_mask)
    [psl, idx] = max(pks(s_mask));
    s_locs = locs(s_mask);
    plot(s_locs(idx), psl, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    text(s_locs(idx)+0.05, psl+2, sprintf('PSL_{v}=%.1f dB', psl), 'Color', 'r', 'FontWeight', 'bold');
end