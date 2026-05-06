%% 稀疏平面阵列设计convex optimazation
% =============================================================
%改10版，无法保存工作空间，最终阵元数为462
% =============================================================
%研二下改11版，第二个凸优化取消阵元坐标
%程序不好？仿照IDEA写？
% =============================================================
clear; clc; close all;
if exist('cvx_setup.m', 'file'), cvx_clear; end

%% --- 1. 参数初始化 ---
%阵列相关
N_side = 32;               % 32x32 候选阵元
N_total = N_side^2;
lamda = 1;                 
d_min = 0.35 * lamda;       % 最小间距约束
k = 2 * pi / lamda;        
Space=lamda/2;             % 阵元初始间距
%性能相关
u1=0.1;                   % 主瓣宽度

rho_mainlobe=0.9;                                %主瓣约束
rho_sll_limit_dB=-30;%规定主副比
%rho_sll_delta_dB=-34;%增加主副比
rho_sll_limit = 10^(rho_sll_limit_dB/20)*rho_mainlobe;        % 旁瓣上限 rho ,对应旁瓣不高于-20dB


delta_lock=0.005;                                % 锁定点阵元激励
test_u = 0.2; test_v = 0.2;%测试坐标
%其他参数
epsilson=1e-6;
sigma=0.01;
threshold = 1e-8;           % 判定阵元是否存活的阈值
mu=1e-4;                    %加权中的极小值常数
eta = 0.05;        % 增大剔除阈值，加快稀疏速度

%---2. 阵列设置 ---
% 初始位置（均匀平面阵）
x_vec = (-15:16) * Space; % 长度为 32，且包含 0,确保（0,0）处恰好有一个阵元
y_vec = (-15:16) * Space;
[X, Y] = meshgrid(x_vec, y_vec);
%[X, Y] = meshgrid(-(N_side-1)/2*Space:Space: (N_side-1)/2*Space, ...
%                  -(N_side-1)/2*Space:Space: (N_side-1)/2*Space );
pos_x = X(:); pos_y = Y(:);
N_total=length(pos_x);     % 更新长度

% ---3. 采样点 (u, v) 空间设置 ---
spacefrequency=50;
u_vec = linspace(-1, 1, spacefrequency); %角域1度的分辨率等效到u/v域，则u/v域约为115个网格点
v_vec = linspace(-1, 1, spacefrequency);
[U, V] = meshgrid(u_vec, v_vec);
u_flat = U(:); v_flat = V(:);
% 法向应指向（0,0），对应的（u，v）也是（0,0）
% 旁瓣区域，32*32的阵列半功率波束宽度在0.05~0.1之间，设为0.12留出裕量
UV_dist=sqrt(u_flat.^2+v_flat.^2);
is_sll = abs(UV_dist) >= u1;
%主瓣区域
is_main = ~is_sll;
center_idx=find(u_flat==0&v_flat==0);%主瓣最高点

%理想主瓣相关设置
delta_m = 0.05; 
alpha = 0.8;
sigma_u = u1/2;
% %主瓣区域和旁瓣区域可视化
% figure('Name', 'UV Region Visualization');
% hold on;
% % 1. 绘制旁瓣采样点 (蓝色)
% scatter(u_flat(is_sll), v_flat(is_sll), 15, [0.8 0.9 1], 'filled', 'DisplayName', 'Sidelobe Region');
% % 2. 绘制主瓣采样点 (红色)
% scatter(u_flat(is_main), v_flat(is_main), 25, 'r', 'filled', 'DisplayName', 'Mainlobe Region');
% % 4. 装饰
% xlabel('u = sin\theta cos\phi');
% ylabel('v = sin\theta sin\phi');
% title(['UV 区域划分 (u_1 = ', num2str(u1), ')']);
% legend('Location', 'northeastoutside');
% axis equal; grid on;
% xlim([-1.1 1.1]); ylim([-1.1 1.1]);
% hold off;
% % --- 打印统计 ---
% fprintf('总采样点数: %d\n', length(u_flat));
% fprintf('主瓣点数: %d\n', sum(is_main));
% fprintf('旁瓣点数: %d\n', sum(is_sll)); 

% %下采样设置
% sll_all_idx = find(is_sll);
% sample_step = 3; 
% sampled_sll_idx = sll_all_idx(1:sample_step:end);%稀疏了1/3

%---4. 迭代前初始化 ---
max_iter = 100;             % 提高最大迭代次数，靠收敛条件触发退出
% 信任域初始化 
delta_w = 0.05;             % 激励信任域半径 δw 
delta_pos =  lamda/60;      % 坐标信任域半径 Δ(k) 
%切比雪夫展开
A_limit = k* delta_pos; % 归一化区间上限
J0=besselj(0, A_limit);
J1=besselj(1, A_limit);
c0 = J0; 
c1 = (2j * J1) ;

%变量设置
w = zeros(N_total, 1);      % 初始激励向量 
num_elements_history = [];  % 记录历史用于绘图
stop_counter = 0;           % 稳定计数器
delta_AF=[];

%扰动设置
gamma=0.05;%扰动幅度比
pos_noise_scale = 0.01 * lamda;
%圆形阵列初始化
R_in = ((N_side-1) * Space) / 2; %圆形半径
offset = Space / 2;
dist_sq = (pos_x - offset).^2 + (pos_y - offset).^2;%由于采用了不对称分布，故将新的圆心加入
in_circle_mask = sqrt(dist_sq) <= (R_in + 1e-6);
w(in_circle_mask) = 1;%圆形内的阵元初始激励设置为1
fprintf('总网格点数: %d, 开启的圆形阵元数: %d\n', N_total, sum(in_circle_mask));
%显示圆形阵列
% figure;
% scatter(pos_x, pos_y, 25, abs(w) ,'filled'); hold on;

% 保证阵元位置不超过矩形
x_min_limit = min(pos_x);
x_max_limit = max(pos_x);
y_min_limit = min(pos_y);
y_max_limit = max(pos_y);

%% --- 2. 迭代优化循环 ---
for iter=1:max_iter
    %迭代开始的更新设置    
    active_count_old =N_total ; % 记录上一轮的阵元数

    %重新计算位置
    dist_from_center=sqrt(pos_x.^2 + pos_y.^2);
    % 锁定矩形阵列和坐标轴交点阵元（不更新位置）
    % 找到边界坐标的极值
    x_max = max(pos_x); x_min = min(pos_x);
    y_max = max(pos_y); y_min = min(pos_y);
    % 1. X轴正半轴交点 (y=0, x=max)
    [~, idx_right] = min(abs(pos_y) + abs(pos_x - x_max));
    % 2. X轴负半轴交点 (y=0, x=min)
    [~, idx_left]  = min(abs(pos_y) + abs(pos_x - x_min));
    % 3. Y轴正半轴交点 (x=0, y=max)
    [~, idx_top]   = min(abs(pos_x) + abs(pos_y - y_max));
    % 4. Y轴负半轴交点 (x=0, y=min)
    [~, idx_bottom]= min(abs(pos_x) + abs(pos_y - y_min));
    % 合并锁定点索引
    lock_idx = [idx_right; idx_left; idx_top; idx_bottom];

    %阵列中心点重新计算
    [~, main_idx] = min(dist_from_center); % 找到初始最接近原点的点

    %加入阵元激励以及坐标的扰动
    if mod(iter, 6) == 0
        fprintf('>>> 触发随机扰动：正在为激励和坐标注入抖动...\n');
        % 1. 激励扰动 (对幅度进行 ±5% 左右的随机抖动，并注入微小相位)
        w_noise_amp = gamma * max(abs(w)); 
        w = w + (randn(size(w)) + 1j*randn(size(w))) * (w_noise_amp / 2);
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
    %% --- 步骤 A: 优化激励 w (固定位置) ---
    AF_base_phase = exp(1j * k * (pos_x * u_flat.' + pos_y * v_flat.'));%单位为1024*13225
    AF0_A=exp(1j * k * (pos_x *0+ pos_y *0));
    Gi = 1 ./ (abs(w) + mu);    %加权,此处的权重应该是上一轮的
    A_ref = exp(-(UV_dist(is_main).^2)/(2*sigma_u^2));
    %检查锁定点与主瓣约束是否冲突
    if any(ismember(main_idx, lock_idx))
    error('主瓣约束点与锁定点发生冲突');
    end
  
    cvx_clear;
    cvx_begin quiet
    variable w_new(N_total, 1) complex
    minimize(norm(Gi.*w_new, 1)) 
    AF_A=AF_base_phase.' * w_new;
    subject to
%         %约束：所有主瓣点的实部为1，虚部为0 (或者根据需求设定范围)（迭代失败，阵元数目没有减小/求解后阵元数目太少）
%         real(AF_A(center_idx))>=1;
%         imag(AF_A(center_idx))==0;

%         %像IDEA一样严格写成u=0&v=0
%         real(AF0_A'*w_new)>=1;
%         imag(AF0_A'*w_new)==0;

        %期望主瓣模版约束√√√√
        real(AF_A(is_main)) <= A_ref + delta_m;
        real(AF_A(is_main)) >= A_ref - delta_m;
        % 限制虚部，确保相位对齐，防止模值因为虚部抵消而凹陷
        imag(AF_A(is_main)) <= delta_m;
        imag(AF_A(is_main)) >= -delta_m;

        % --- 2. 信任域约束 ---
        % norm(w_new - w, 2) <= delta_w; 
        % --- 3. 旁瓣约束（下采样） ---
        %abs(AF_A(is_sll))<= rho_sll_limit;
        norm(AF_A(is_sll),inf)<=rho_sll_limit;
        %锁定点约束
        w_new(lock_idx)==delta_lock;

cvx_end
    % 如果 CVX 失败，跳出循环
    if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
        fprintf('\n[警告] 迭代%d阵元激励求解失败。',iter);
        break;
    else
        keep_idx = abs(w_new) >  eta * max(abs(w_new));
        if sum(keep_idx)<7
            fprintf('\n[警告]迭代%d阵元激励求解后阵元太少，停止迭代。',iter);
            break;
        else
        w = w_new;%更新变量
        fprintf('Iteration %d:阵元激励优化Status = %s\n', iter,  cvx_status);
        end
    end
            
%     %% --- 步骤 B: 优化坐标增量 (固定激励 w_new, 优化位置以压低旁瓣) ---
%     %优化目标改为最小化主瓣区域与理想高斯响应的区别
%     cvx_clear;
%     cvx_begin quiet
%         variables dx(N_total, 1) dy(N_total, 1)    
%         xi_val = k * (dx * u_flat.' + dy * v_flat.');
%         approx_exp = c0 + (c1/ A_limit) * xi_val ;
%         AF_B = (w.' * (AF_base_phase .* approx_exp)).';
%         %对于主瓣
%         approx_B=c0 + (c1/ A_limit) *(k * (dx * 0 + dy * 0));
%         AF0_B=(w.' * (AF0_A .* approx_B)).';
%         % 目标函数
%         minimize(norm(real(AF_B(is_main)) - A_ref, 2))
%         subject to
% %             % 约束更新后的位置不越界
% %             %此处的最大最小限制为循环外定义，所以是最初的矩形阵列边界，不会随着循环而改变
% %             pos_x + dx >= x_min_limit;
% %             pos_x + dx <= x_max_limit;
% %             pos_y + dy >= y_min_limit;
% %             pos_y + dy <= y_max_limit;
%             %锁定点位置不改变
%             dx(lock_idx) == 0; 
%             dy(lock_idx) == 0; 
%             %信任域约束
%             norm([dx; dy], 2) <= delta_pos; % 坐标信任域约束
%            
%             % 旁瓣约束：
%             %abs(AF_B(is_sll))<= t_sll;
%             norm(AF_B(is_sll), inf) <=rho_sll_limit ;
% 
%           %主瓣约束
%           real(AF_B(center_idx)) >= rho_mainlobe; 
%           imag(AF_B(center_idx)) == 0;
% 
% %             %像IDEA一样写成u=0&v=0
% %             real(AF0_B)>=rho_mainlobe
% %             imag(AF0_B)==0;
% 
%             % 4. 辅助约束：限制主瓣虚部，防止相位剧烈旋转导致幅度抵消
%             abs(imag(AF_B(is_main))) <= delta_m;
%              %期望主瓣模版约束
% %             real(AF_B(is_main)) <= A_ref + delta_m;
% %             real(AF_B(is_main)) >= A_ref - delta_m;
% %             % 限制虚部，确保相位对齐，防止模值因为虚部抵消而凹陷
% %             imag(AF_B(is_main)) <= delta_m;
% %             imag(AF_B(is_main)) >= -delta_m;
%             
%                 % 阵元间距约束 (线性化处理：保证不重叠)
%                 for i = 1:N_total
%                 % 只计算 i 之后的阵元，避免重复约束
%                 for j = i+1:N_total
%                     % 当前实际物理距离
%                     dist_x = pos_x(i) - pos_x(j);
%                     dist_y = pos_y(i) - pos_y(j);
%                     current_dist = sqrt(dist_x^2 + dist_y^2);
%                     
%                     % 如果两阵元靠得太近，则施加线性化排斥约束
%                     if current_dist <  d_min
%                         % 对非凸的阵元间距约束进行泰勒一阶展开
%                         % 约束形式：current_dist + (dist_x*(dx(i)-dx(j)) + dist_y*(dy(i)-dy(j)))/current_dist >= d_min
%                         (dist_x * (dx(i) - dx(j)) + dist_y * (dy(i) - dy(j))) / current_dist+epsilson >= d_min - current_dist;
%                     end
%                 end
%             end
% %                % --- 优化后的间距约束处理 ---
% %                 for i = 1:N_total
% %                     for j = i+1:N_total
% %                         dist_x = pos_x(i) - pos_x(j);
% %                         dist_y = pos_y(i) - pos_y(j);
% %                         current_dist = sqrt(dist_x^2 + dist_y^2);
% %                         
% %                         % 只有当距离小于最小间距且大于一个极小安全距离时才添加约束
% %                         if current_dist < d_min && current_dist > 1e-4
% %                             % 将除法移项或者确保分母安全
% %                             % 线性化公式: (dist_x*ddx + dist_y*ddy) / dist >= d_min - dist
% %                             inv_dist = 1.0 / current_dist;
% %                             (dist_x * (dx(i) - dx(j)) + dist_y * (dy(i) - dy(j))) * inv_dist >= d_min - current_dist;
% %                         elseif current_dist <= 1e-4
% %                             % 如果阵元完全重叠，手动强制它们在本次迭代向相反方向移动
% %                             dx(i) - dx(j) >= 0.01 * d_min;
% %                             dy(i) - dy(j) >= 0.01 * d_min;
% %                         end
% %                     end
% %                 end
%     cvx_end
%     
%     % 容错处理
%     %如果 CVX 失败，跳出循环
%     if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
%         fprintf('警告：第%d次迭代坐标优化失败，保持位置不变。\n',iter);
%         dx = zeros(N_total, 1); 
%         dy = zeros(N_total, 1); 
%         break;
%     else
%         %更新变量
%         pos_x=pos_x+dx;
%         pos_y=pos_y+dy;
%         fprintf('Iteration %d:阵元坐标优化Status = %s\n', iter,  cvx_status);
%     end
    %% --- 步骤 B: 基于扩展正多边形采样与激励重选的坐标更新 ---
% M_poly = 6;           % 扩展正多边形边数（如六边形）
% R_search = delta_pos; % 搜索半径，复用原信任域半径
% 
% % 1. 生成候选坐标矩阵 (N_total x (M+1))
% % 每行包含：[原坐标, 顶点1, 顶点2, ..., 顶点M]
% candidate_x = zeros(N_total, M_poly + 1);
% candidate_y = zeros(N_total, M_poly + 1);
% 
% for i = 1:N_total
%     candidate_x(i, 1) = pos_x(i);
%     candidate_y(i, 1) = pos_y(i);
%     for m = 1:M_poly
%         theta = 2 * pi * (m-1) / M_poly;
%         candidate_x(i, m+1) = pos_x(i) + R_search * cos(theta);
%         candidate_y(i, m+1) = pos_y(i) + R_search * sin(theta);
%     end
% end
% 
% % 展平候选点以便计算流形矩阵 (Total_Candidates = N_total * (M+1))
% flat_cand_x = candidate_x(:);
% flat_cand_y = candidate_y(:);
% N_cand_total = length(flat_cand_x);
% 
% % 2. 重新构建流形矩阵 (针对所有候选位置)
% AF_cand_phase = exp(1j * k * (flat_cand_x * u_flat.' + flat_cand_y * v_flat.'));
% 
% cvx_clear;
% cvx_begin quiet
%     % 变量为所有候选点的激励
%     variable w_cand(N_cand_total, 1) complex
%     
%     % 构造当前位置的等效方向图
%     AF_B = (AF_cand_phase.' * w_cand); 
%     
%     % 优化目标：最小化主瓣误差 + 稀疏项
%     % 注意：此处引入 Gi_cand 可以进一步引导稀疏
%     minimize( norm(real(AF_B(is_main)) - A_ref, 2) + 0.1 * norm(w_cand, 1) )
%     
%     subject to
%         % 旁瓣约束
%         norm(AF_B(is_sll), inf) <= rho_sll_limit;
%         
%         % 主瓣中心约束
%         real(AF_B(center_idx)) >= rho_mainlobe;
%         imag(AF_B(center_idx)) == 0;
%         
%         % 间距约束 (近似处理：只对原位置进行锁定点保护)
%         % 实际中候选点过多，全局间距约束会导致计算极慢，建议仅对锁定点强制约束
%         % 或者在更新坐标后进行物理碰撞检查
%         w_cand((lock_idx-1)*(M_poly+1) + 1) == delta_lock; 
%         
%         % (可选) 强制每个阵元的候选组中只有一个能存活，
%         % 但为了计算效率，通常靠 minimize(norm(w_cand,1)) 自动竞争
% cvx_end
% 
% if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
%     fprintf('警告：第%d次迭代坐标采样优化失败。\n', iter);
% else
%     % 3. 坐标重选逻辑
%     new_pos_x = zeros(N_total, 1);
%     new_pos_y = zeros(N_total, 1);
%     new_w = zeros(N_total, 1);
%     
%     w_reshaped = reshape(w_cand, M_poly + 1, N_total).'; % 变为 N_total x (M+1)
%     
%     for i = 1:N_total
%         % 寻找该阵元所有候选点中激励模值最大的索引
%         [max_val, best_cand_idx] = max(abs(w_reshaped(i, :)));
%         
%         % 更新为最佳候选点的坐标和激励
%         new_pos_x(i) = candidate_x(i, best_cand_idx);
%         new_pos_y(i) = candidate_y(i, best_cand_idx);
%         new_w(i) = w_reshaped(i, best_cand_idx);
%     end
%     
%     pos_x = new_pos_x;
%     pos_y = new_pos_y;
%     w = new_w;
%     fprintf('Iteration %d: 坐标多边形重选完成，Status = %s\n', iter, cvx_status);
% end
%     %% --- 步骤 C: 计算线性近似损失  ---
%     % 选取测试点评估 Pade 近似预测效果
%     % Pade 近似预测值
%     xi_test =  k * (dx * test_u + dy * test_v);
%     P_plus = sum(w .* exp(1j*k*(pos_x*test_u + pos_y*test_v)) .* (c0 + (c1/ A_limit) * xi_test));
%     % 实际阵列因子值 
%     AF_plus = sum(w .* exp(1j*k*((pos_x+dx)*test_u + (pos_y+dy)*test_v)));
%     %使用近似之后的方向图和原方向图的差值
%     delta_AF(end+1)=P_plus-AF_plus;
%   
%% 停止迭代以及更新
    active_idx = abs(w) > threshold;
    active_idx(lock_idx) = true;       % 强制将锁定点设为 true 
    active_count_new = sum(active_idx);
    if active_count_new<7
        fprintf('\n[警告] 阵元太少，停止本轮迭代。');
        break; % 阵元太少则停止
    else
        num_elements_history(end+1) = active_count_new;
        pos_x=pos_x(active_idx);
        pos_y=pos_y(active_idx);
        w=w(active_idx);
        fprintf('Iteration %d: Active Elements = %d\n', iter, active_count_new);
    end
    %阵元激励向量长度更新
    N_total=length(w);
    
    %停止迭代判定逻辑
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
%% 函数引用部分
%显示坐标换成1000
Axis_scale=1000;
%满阵加窗比较
[U_dB_taylor,U_dB_hamming,U_dB_hanning,U_dB_ChebyShev,V_dB_taylor,V_dB_hamming,V_dB_hanning,V_dB_ChebyShev...
    ] = full_array_with_windows(N_side,k,Space,rho_sll_limit_dB,Axis_scale);
%阵元激励为全1的满阵比较
 [U_dB_full_array,V_dB_full_array] = full_array(N_side,k,Space,Axis_scale);

% %% 可视化线性近似随迭代的损失
% figure;plot(1:length(delta_AF),delta_AF,'-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
% title('线性近似的损失随迭代变化过程'); xlabel('迭代次数'); ylabel('线性近似的损失'); grid on;

%% --- 3. 可视化 ---

% 1. 阵元数目随迭代的变化
figure; plot(num_elements_history, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
title('阵元数目随迭代变化过程'); xlabel('迭代次数'); ylabel('活跃阵元数'); grid on;

% 2. 最终阵元分布
figure; 
scatter(pos_x, pos_y, 25, abs(w) ,'filled'); hold on;
[~, idx_r] = min(abs(pos_y) + abs(pos_x - x_max_limit));
[~, idx_l] = min(abs(pos_y) + abs(pos_x - x_min_limit)); 
[~, idx_t] = min(abs(pos_x) + abs(pos_y - y_max_limit)); 
[~, idx_b] = min(abs(pos_x) + abs(pos_y - y_min_limit)); 
lock_original_idx=[idx_r,idx_l,idx_t,idx_b];
scatter(pos_x(lock_original_idx), pos_y(lock_original_idx), 100, 'r', 'x', 'LineWidth', 2);
title('优化后的阵元分布 (红叉为锁定角点)'); axis equal; xlabel('x(m)'); ylabel('y(m)');

% 3. 三维方向图
slice_idx_u= linspace(-1,1,Axis_scale);
slice_idx_v= linspace(-1,1,Axis_scale);
[Slice_U,Slice_V]=meshgrid(slice_idx_u,slice_idx_v);
slice_u=Slice_U(:);slice_v=Slice_V(:);
Phase_Mat = exp(1j * k * (pos_x* slice_u.' + pos_y * slice_v.'));
AF_show = w.' * Phase_Mat;
[rows,cols]=size(Slice_U);
AF_amp  = reshape(AF_show,rows,cols);
AF_dB = 20 * log10(abs(AF_amp)/max(abs(AF_amp(:))));%转dB，做归一化

% --- 4. 方向图绘图 ---
figure;
surf(slice_idx_u, slice_idx_v, AF_dB);
shading interp;
colormap(jet);
colorbar;
view(35, 35);
axis tight;
xlabel('u');
ylabel('v');
zlabel('幅度 (dB)');
title('优化后三维方向图（u,v）');

%% --- u 方向切片（v = 0） ---
AF_u_slice =w.' * exp(1j * k * (pos_x* slice_idx_u));
% 预先转换为 dB 量纲，方便寻峰与标注
AF_u_dB = 20 * log10(abs(AF_u_slice)/max(abs(AF_u_slice(:))) ); 
figure;
plot(slice_idx_u, AF_u_dB, 'LineWidth', 1.8,'DisplayName', '稀疏优化阵列');
hold on; 
plot(slice_idx_u, U_dB_taylor, 'r--', 'LineWidth', 1.2, 'DisplayName', '满阵: Taylor Window');
plot(slice_idx_u, U_dB_ChebyShev, 'g-.', 'LineWidth', 1.2, 'DisplayName', '满阵: Chebyshev Window');
plot(slice_idx_u, U_dB_hamming, 'k:', 'LineWidth', 1.2, 'DisplayName', '满阵: Hamming Window');
plot(slice_idx_u, U_dB_hanning, 'm:', 'LineWidth', 1.2, 'DisplayName', '满阵: Hanning Window');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on;
ylim([-40 5]); 
xlim([min(slice_idx_u) max(slice_idx_u)]);
xlabel('u = sin(\theta)cos(\phi)');
ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）');
% 4. 旁瓣标注逻辑 (基于 dB 数据寻峰)
[all_pks, all_locs] = findpeaks(AF_u_dB, slice_idx_u);
% 识别并剔除主瓣 (通常是最大的峰值)
[~, mainpeak_idx] = max(all_pks); 
if ~isempty(mainpeak_idx)
    all_pks(mainpeak_idx) = -inf; 
end
% 寻找最高旁瓣
[max_sll_val, max_sll_idx] = max(all_pks);
u_max_sll = all_locs(max_sll_idx);
% 5. 绘制标注
if ~isempty(max_sll_idx)
    % 绘制红圈标记
    plot(u_max_sll, max_sll_val, 'ro', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility', 'off');
    % 动态调整文字位置
    if u_max_sll > 0
        align = 'left'; offset = 0.05; 
    else
        align = 'right'; offset = -0.05; 
    end
    % 显示标注文字
    text(u_max_sll + offset, max_sll_val, ...
        sprintf('  最高旁瓣: %.2f dB\n  u = %.4f', max_sll_val, u_max_sll), ...
        'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', align);
    % 绘制一条水平虚线显示旁瓣水平 (可选，增加可读性)
    line([min(slice_idx_u) max(slice_idx_u)], [max_sll_val max_sll_val], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1,'HandleVisibility', 'off');
end
hold off;

%分别计算3dB波束宽度
[hpbw30_ChebyshevSCP_u] = calculate_hpbw(slice_idx_u, AF_u_dB);
[hpbw30_taylor_u] = calculate_hpbw(slice_idx_u, U_dB_taylor);
[hpbw30_ChebyShev_u] = calculate_hpbw(slice_idx_u,U_dB_ChebyShev);
[hpbw30_hamming_u] = calculate_hpbw(slice_idx_u,U_dB_hamming);
[hpbw30_hanning_u] = calculate_hpbw(slice_idx_u,U_dB_hanning);

fprintf("--------------------方向维-30dB主瓣宽度比较--------------------\n");
fprintf('u轴稀疏阵列 (3dB宽度): %.4f°\n', hpbw30_ChebyshevSCP_u);
fprintf('u轴满阵加泰勒窗 (3dB宽度): %.4f°\n', hpbw30_taylor_u);
fprintf('u轴满阵加切比雪夫窗 (3dB宽度): %.4f°\n', hpbw30_ChebyShev_u);
fprintf('u轴满阵加汉明窗 (3dB宽度): %.4f°\n', hpbw30_hamming_u);
fprintf('u轴满阵加汉宁窗 (3dB宽度): %.4f°\n', hpbw30_hanning_u);

%% % ==============u方向切片和满阵比较主瓣宽度和主副比======================================
figure;
plot(slice_idx_u, AF_u_dB, 'LineWidth', 1.8,'DisplayName', '稀疏优化阵列');
hold on; 
plot(slice_idx_u, U_dB_full_array, 'm--', 'LineWidth', 1.8, 'DisplayName', '满阵（阵元激励为1）');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴基本属性
grid on;
ylim([-40 5]); 
xlim([min(slice_idx_u) max(slice_idx_u)]);
xlabel('u = sin(\theta)cos(\phi)');
ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）');

% 4.稀疏阵列优化旁瓣标注逻辑 (基于 dB 数据寻峰)
if ~isempty(max_sll_idx)% 绘制标注
    plot(u_max_sll, max_sll_val, 'ro', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility', 'off'); % 绘制红圈标记
    if u_max_sll > 0% 动态调整文字位置
        align = 'left'; offset = 0.3; 
    else
        align = 'right'; offset = -0.3; 
    end
    text(u_max_sll + offset, max_sll_val, ...
        sprintf('  最高旁瓣: %.2f dB\n  u = %.4f', max_sll_val, u_max_sll), ...
        'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', align);% 显示标注文字
    line([min(slice_idx_u) max(slice_idx_u)], [max_sll_val max_sll_val], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1,'HandleVisibility', 'off');% 绘制一条水平虚线显示旁瓣水平 (可选，增加可读性)
end
% 5. 满阵旁瓣标注逻辑 (基于 dB 数据寻峰)
[all_pks, all_locs] = findpeaks(U_dB_full_array, slice_idx_u);
[~, mainpeak_idx] = max(all_pks); % 识别并剔除主瓣 (通常是最大的峰值)
if ~isempty(mainpeak_idx)
    all_pks(mainpeak_idx) = -inf; 
end
[max_sll_val2, max_sll_idx2] = max(all_pks);% 寻找最高旁瓣
u_max_sll2 = all_locs(max_sll_idx2);
if ~isempty(max_sll_idx2)% 5. 绘制标注
    plot(u_max_sll2, max_sll_val2, 'ro', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility', 'off'); % 绘制红圈标记
    if u_max_sll2 > 0% 动态调整文字位置
        align = 'left'; offset = 0.05; 
    else
        align = 'right'; offset = -0.05; 
    end
    text(u_max_sll2 + offset, max_sll_val2, ...
        sprintf('  最高旁瓣: %.2f dB\n  u = %.4f', max_sll_val2, u_max_sll2), ...
        'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', align); % 显示标注文字
    line([min(slice_idx_u) max(slice_idx_u)], [max_sll_val2 max_sll_val2], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1,'HandleVisibility', 'off'); % 绘制一条水平虚线显示旁瓣水平 (可选，增加可读性)
end
hold off;

%% --- v 方向切片（u = 0） ---
AF_v_slice=w.' * exp(1j * k * (pos_y* slice_idx_v));
% 2. 转换为 dB 量纲
AF_v_dB = 20 * log10(abs(AF_v_slice)/max(abs(AF_v_slice(:))));
figure;
% 使用橙红色系曲线以区分 u 轴切片
plot(slice_idx_v, AF_v_dB, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980],'DisplayName', '稀疏优化阵列'); 
hold on;
plot(slice_idx_v, V_dB_taylor, 'r--', 'LineWidth', 1.2, 'DisplayName', '满阵: Taylor Window');
plot(slice_idx_v, V_dB_ChebyShev, 'g-.', 'LineWidth', 1.2, 'DisplayName', '满阵: Chebyshev Window');
plot(slice_idx_v, V_dB_hamming, 'k:', 'LineWidth', 1.2, 'DisplayName', '满阵: Hamming Window');
plot(slice_idx_v, V_dB_hanning, 'm:', 'LineWidth', 1.2, 'DisplayName', '满阵: Hanning Window');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴属性
grid on;
ylim([-40 5]); 
xlim([min(slice_idx_v) max(slice_idx_v)]);
xlabel('v = sin(\theta)sin(\phi)');
ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)');
% 4. 旁瓣标注逻辑 (基于 dB 数据寻峰)
[all_pks_v, all_locs_v] = findpeaks(AF_v_dB, slice_idx_v);
% 识别并剔除主瓣 (峰值最大处)
[~, mainpeak_idx_v] = max(all_pks_v); 
if ~isempty(mainpeak_idx_v)
    all_pks_v(mainpeak_idx_v) = -inf; 
end
% 寻找最高旁瓣
[max_v_sll_val, max_v_sll_idx] = max(all_pks_v);
v_max_sll_pos = all_locs_v(max_v_sll_idx);
% 5. 绘制标注
if ~isempty(max_v_sll_idx)
    % 绘制红圈标记最高旁瓣点
    plot(v_max_sll_pos, max_v_sll_val, 'ro', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility', 'off');
    % 根据位置动态调整文字对齐方式
    if v_max_sll_pos > 0
        v_align = 'left'; v_offset = 0.05; 
    else
        v_align = 'right'; v_offset = -0.05; 
    end
    % 显示标注文字 (显示 dB 值和对应的 v 坐标)
    text(v_max_sll_pos + v_offset, max_v_sll_val, ...
        sprintf('  v轴最高旁瓣: %.2f dB\n  v = %.4f', max_v_sll_val, v_max_sll_pos), ...
        'Color', [0.5 0 0], 'FontWeight', 'bold', 'HorizontalAlignment', v_align);
    % 绘制水平参考线
    line([min(slice_idx_v) max(slice_idx_v)], [max_v_sll_val max_v_sll_val], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1,'HandleVisibility', 'off');
end
hold off;

%分别计算3dB波束宽度
[hpbw30_ChebyshevSCP_v] = calculate_hpbw(slice_idx_v, AF_v_dB);
[hpbw30_taylor_v] = calculate_hpbw(slice_idx_v, V_dB_taylor);
[hpbw30_ChebyShev_v] = calculate_hpbw(slice_idx_v,V_dB_ChebyShev);
[hpbw30_hamming_v] = calculate_hpbw(slice_idx_v,V_dB_hamming);
[hpbw30_hanning_v] = calculate_hpbw(slice_idx_v,V_dB_hanning);

fprintf("--------------------俯仰维-30dB主瓣宽度比较--------------------\n");
fprintf('v轴稀疏阵列 (3dB宽度): %.4f°\n', hpbw30_ChebyshevSCP_v);
fprintf('v轴满阵加泰勒窗 (3dB宽度): %.4f°\n', hpbw30_taylor_v);
fprintf('v轴满阵加切比雪夫窗 (3dB宽度): %.4f°\n', hpbw30_ChebyShev_v);
fprintf('v轴满阵加汉明窗 (3dB宽度): %.4f°\n', hpbw30_hamming_v);
fprintf('v轴满阵加汉宁窗 (3dB宽度): %.4f°\n', hpbw30_hanning_v);

%% ==========v方向切片和满阵比较主瓣宽度和主副比===============================
figure;
% 使用橙红色系曲线以区分 u 轴切片
plot(slice_idx_v, AF_v_dB, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980],'DisplayName', '稀疏优化阵列'); 
hold on;
plot(slice_idx_v, V_dB_full_array, 'b--', 'LineWidth', 1.8, 'DisplayName', '满阵（阵元激励为1）');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴属性
grid on;
ylim([-40 5]); 
xlim([min(slice_idx_v) max(slice_idx_v)]);
xlabel('v = sin(\theta)sin(\phi)');
ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)');
% 4. 稀疏优化阵列旁瓣标注逻辑 (基于 dB 数据寻峰)
if ~isempty(max_v_sll_idx)%  绘制标注
    plot(v_max_sll_pos, max_v_sll_val, 'ro', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility', 'off');% 绘制红圈标记最高旁瓣点
    if v_max_sll_pos > 0 % 根据位置动态调整文字对齐方式
        v_align = 'left'; v_offset = 0.3; 
    else
        v_align = 'right'; v_offset = -0.3; 
    end
    text(v_max_sll_pos + v_offset, max_v_sll_val, ...
        sprintf('  v轴最高旁瓣: %.2f dB\n  v = %.4f', max_v_sll_val, v_max_sll_pos), ...
        'Color', [0.5 0 0], 'FontWeight', 'bold', 'HorizontalAlignment', v_align);% 显示标注文字 (显示 dB 值和对应的 v 坐标)
    line([min(slice_idx_v) max(slice_idx_v)], [max_v_sll_val max_v_sll_val], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1,'HandleVisibility', 'off');% 绘制水平参考线
end
% 5. 满阵旁瓣标注逻辑 (基于 dB 数据寻峰)
[all_pks_v, all_locs_v] = findpeaks(V_dB_full_array, slice_idx_v);
[~, mainpeak_idx_v] = max(all_pks_v); % 识别并剔除主瓣 (峰值最大处)
if ~isempty(mainpeak_idx_v)
    all_pks_v(mainpeak_idx_v) = -inf; 
end
[max_v_sll_val2, max_v_sll_idx2] = max(all_pks_v);% 寻找最高旁瓣
v_max_sll_pos2 = all_locs_v(max_v_sll_idx2);
if ~isempty(max_v_sll_idx2)%  绘制标注
    plot(v_max_sll_pos2, max_v_sll_val2, 'ro', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility', 'off');% 绘制红圈标记最高旁瓣点
    if v_max_sll_pos2 > 0 % 根据位置动态调整文字对齐方式
        v_align = 'left'; v_offset = 0.05; 
    else
        v_align = 'right'; v_offset = -0.05; 
    end
    text(v_max_sll_pos2 + v_offset, max_v_sll_val2, ...
        sprintf('  v轴最高旁瓣: %.2f dB\n  v = %.4f', max_v_sll_val2, v_max_sll_pos2), ...
        'Color', [0.5 0 0], 'FontWeight', 'bold', 'HorizontalAlignment', v_align);% 显示标注文字 (显示 dB 值和对应的 v 坐标)
    line([min(slice_idx_v) max(slice_idx_v)], [max_v_sll_val2 max_v_sll_val2], ...
        'Color', 'r', 'LineStyle', '--', 'LineWidth', 1,'HandleVisibility', 'off');% 绘制水平参考线
end
hold off;

%% 相同副瓣的满阵加窗比较
 %u方向函数引用
[Slice_U_taylor,Slice_U_ChebyShev] = full_array_with_windows_aftersamesidelobe_u(N_side,k,Space,max_sll_val,Axis_scale);
 %v方向函数引用
[Slice_V_taylor,Slice_V_ChebyShev] = full_array_with_windows_aftersamesidelobe_v(N_side,k,Space,max_v_sll_val,Axis_scale);

%u轴切片方向图显示
figure;
plot(slice_idx_u, AF_u_dB, 'LineWidth', 1.8,'DisplayName', '稀疏优化阵列');
hold on; % 
plot(slice_idx_u, Slice_U_taylor, 'm--', 'LineWidth', 1.8, 'DisplayName', '满阵加泰勒窗');
plot(slice_idx_u, Slice_U_ChebyShev, 'c--', 'LineWidth', 1.8, 'DisplayName', '满阵加切比雪夫窗');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴基本属性
grid on;
ylim([-40 5]); % 建议上限给到 5dB，防止主瓣顶端被压住
xlim([min(slice_idx_u) max(slice_idx_u)]);
xlabel('u = sin(\theta)cos(\phi)');
ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）（相同旁瓣电平）');
hold off;

%v轴切片方向图显示
figure;
% 使用橙红色系曲线以区分 u 轴切片
plot(slice_idx_v, AF_v_dB, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980],'DisplayName', '稀疏优化阵列'); 
hold on;
plot(slice_idx_v, Slice_V_taylor, 'b--', 'LineWidth', 1.8, 'DisplayName', '满阵加泰勒窗');
plot(slice_idx_v, Slice_V_ChebyShev, 'y--', 'LineWidth', 1.8, 'DisplayName', '满阵加切比雪夫窗');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
% 3. 设置坐标轴属性
grid on;
ylim([-40 5]); 
xlim([min(slice_idx_v) max(slice_idx_v)]);
xlabel('v = sin(\theta)sin(\phi)');
ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)（相同旁瓣电平）');

%分别计算3dB波束宽度
[hpbw_sparse_u] = calculate_hpbw(slice_idx_u, AF_u_dB);
[hpbw_taylor_u] = calculate_hpbw(slice_idx_u, Slice_U_taylor);
[hpbw_ChebyShev_u] = calculate_hpbw(slice_idx_u,Slice_U_ChebyShev);

fprintf("----------------相同副瓣的满阵加窗比较------------------------\n");
fprintf('u轴稀疏阵列 (3dB宽度): %.4f°\n', hpbw_sparse_u);
fprintf('u轴满阵加泰勒窗 (3dB宽度): %.4f°\n', hpbw_taylor_u);
fprintf('u轴满阵加切比雪夫窗 (3dB宽度): %.4f°\n', hpbw_ChebyShev_u);

[hpbw_sparse_v] = calculate_hpbw(slice_idx_v, AF_v_dB);
[hpbw_ChebyShev_v] = calculate_hpbw(slice_idx_v, Slice_V_taylor);
[hpbw_taylor_v] = calculate_hpbw(slice_idx_v,Slice_V_ChebyShev);

fprintf('v轴稀疏阵列 (3dB宽度): %.4f°\n', hpbw_sparse_v);
fprintf('v轴满阵加泰勒窗 (3dB宽度): %.4f°\n', hpbw_taylor_v);
fprintf('v轴满阵加切比雪夫窗 (3dB宽度): %.4f°\n', hpbw_ChebyShev_v);
%% 其他方法的对比
%IDEA
folderpath='C:\Users\gyq\Desktop\sparse_arrays_design\others_method\IDEA';
filename='IDEA5WorkSpace.mat';
fullPath=fullfile(folderpath, filename);

dataStruct=load(fullPath,'U_f','V_f','pos_x','pos_y','AF_dB','AF','c','slice_axis','pattern_u','pattern_v','show_sampling');

fprintf('IDEA方法最终阵元数目为：%d\n',length(dataStruct.pos_x));

%u轴切片方向图显示
figure;
plot(slice_idx_u, AF_u_dB, 'LineWidth', 1.8,'DisplayName', '稀疏优化阵列');
hold on; % ！！！关键点：必须 hold on 才能在同一张图上画标注
plot(slice_idx_u,dataStruct.pattern_u, 'm--', 'LineWidth', 1.8, 'DisplayName', 'IDEA');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on;
ylim([-40 5]); % 建议上限给到 5dB，防止主瓣顶端被压住
xlim([min(slice_idx_u) max(slice_idx_u)]);
xlabel('u = sin(\theta)cos(\phi)');
ylabel('幅度 (dB)');
title('u 轴切片方向图（v=0）（不同优化方法比较）');
hold off;
%v轴切片方向图显示
figure;
plot(slice_idx_v, AF_v_dB, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980],'DisplayName', '稀疏优化阵列'); % 使用橙红色系曲线以区分 u 轴切片
hold on;
plot(slice_idx_v, dataStruct.pattern_v, 'b--', 'LineWidth', 1.8, 'DisplayName', 'IDEA');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'k');
grid on;
ylim([-40 5]); 
xlim([min(slice_idx_v) max(slice_idx_v)]);
xlabel('v = sin(\theta)sin(\phi)');
ylabel('幅度 (dB)');
title('v 轴切片方向图 (u=0)（不同优化方法比较）');
%% 和其他优化方法比较方向性系数
U_power = abs(AF_amp).^2;
U_max = max(U_power(:));
du = 2/(Axis_scale-1);
dv = 2/(Axis_scale-1);
valid_mask = (Slice_U.^2 + Slice_V.^2) <= 1;
P_rad = sum(U_power(valid_mask)) * (du * dv); % du, dv 为采样步长
Directivity_linear = (4 * pi * U_max) / P_rad;
Directivity_dB = 10 * log10(Directivity_linear);
fprintf('本稀疏阵列优化方法的方向性系数为: %.2f dB\n', Directivity_dB);

U_power_IDEA = abs(dataStruct.AF).^2;
U_max_IDEA = max(U_power_IDEA(:));
du_IDEA = 2/(dataStruct.show_sampling-1);
dv_IDEA = 2/(dataStruct.show_sampling-1);
valid_mask_IDEA = (dataStruct.U_f.^2 + dataStruct.V_f.^2) <= 1;
P_rad_IDEA = sum(U_power_IDEA(valid_mask_IDEA)) * (du_IDEA * dv_IDEA); % du, dv 为采样步长
Directivity_linear_IDEA = (4 * pi * U_max_IDEA) / P_rad_IDEA;
Directivity_dB_IDEA = 10 * log10(Directivity_linear_IDEA);
fprintf('IDEA优化方法的方向性系数为: %.2f dB\n', Directivity_dB_IDEA);
