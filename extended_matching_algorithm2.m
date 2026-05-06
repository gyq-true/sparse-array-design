function [dx_recovered, dy_recovered] = extended_matching_algorithm2(old_pos_x, old_pos_y, Q_size, new_pos_x, new_pos_y, AF_B, u_flat, v_flat, w, k)
    %% 改进的唯一性匹配算法 (Unique Pairing with d_min Penalty)
    %GEMINI改，本人未修改
    
    % 参数设置
    lamda = 2 * pi / k;
    d_min = 0.35 * lamda; % 这里的阈值应与主程序一致
    
    % 初始化
    % 确保输入的是列向量
    new_pos_x = new_pos_x(:);
    new_pos_y = new_pos_y(:);
    z_best = new_pos_y; 
    
    % 计算初始误差
    AF_init = (exp(1j * k * (new_pos_x * u_flat.' + z_best * v_flat.'))).' * w;
    current_best_delta = sum(abs(AF_init - AF_B).^2) / sum(abs(AF_B).^2);

    % --- 步骤 1: 迭代交换匹配 (确保 1对1 且误差最小) ---
    max_swaps = 2; % 扫描遍数，防止计算量过大
    for s = 1:max_swaps
        improved = false;
        for m = 1 : (Q_size - 1)
            for n = (m + 1) : Q_size
                % 尝试交换第 m 和第 n 个 y 坐标
                z_temp = z_best;
                z_temp(m) = z_best(n);
                z_temp(n) = z_best(m);
                
                % 计算交换后的阵列因子
                % 优化：只计算受影响的两个阵元的变化（可选），这里为了稳妥采用全量计算
                F_temp = (exp(1j * k * (new_pos_x * u_flat.' + z_temp * v_flat.'))).' * w;
                
                % 计算误差：方向图误差 + 惩罚项（如果太近则误差激增）
                dist_mn = sqrt((new_pos_x(m)-new_pos_x(n))^2 + (z_temp(m)-z_temp(n))^2);
                penalty = 0;
                if dist_mn < d_min
                    penalty = 100; % 给予极大的惩罚，阻止这种匹配
                end
                
                delta_temp = sum(abs(F_temp - AF_B).^2) / sum(abs(AF_B).^2) + penalty;
                
                % 如果误差确实变小且满足约束
                if delta_temp < current_best_delta
                    z_best = z_temp;
                    current_best_delta = delta_temp;
                    improved = true;
                end
            end
        end
        if ~improved, break; end % 如果一整轮都没有改进，提前退出
    end
    
    dx_paired = new_pos_x;
    dy_paired = z_best;

    % --- 步骤 2: 基于匈牙利算法（或唯一性最近邻）恢复原始顺序 ---
    % 这一步保证了优化后的坐标点 (dx_paired, dy_paired) 能够“对号入座”回到旧阵元顺序中
    dx_recovered = zeros(Q_size, 1);
    dy_recovered = zeros(Q_size, 1);
    
    % 构造代价矩阵 (欧氏距离)
    cost_mat = zeros(Q_size, Q_size);
    for i = 1:Q_size
        cost_mat(i, :) = sqrt((old_pos_x(i) - dx_paired).^2 + (old_pos_y(i) - dy_paired).^2);
    end
    
    % 使用 MATLAB 自带的指派函数（需要优化工具箱）
    % 如果没有该工具箱，保留你原有的 used_new_idx 逻辑也是唯一的
    try
        assignment = matchpairs(cost_mat, 100); % 100 为允许的最大代价
        for k_idx = 1:size(assignment, 1)
            dx_recovered(assignment(k_idx,1)) = dx_paired(assignment(k_idx,2));
            dy_recovered(assignment(k_idx,1)) = dy_paired(assignment(k_idx,2));
        end
    catch
        % 如果没有 matchpairs 函数，使用鲁棒的贪婪唯一性匹配
        used_new_idx = false(Q_size, 1);
        for i = 1:Q_size
            current_dist_row = cost_mat(i, :);
            current_dist_row(used_new_idx) = inf;
            [~, best_j] = min(current_dist_row);
            dx_recovered(i) = dx_paired(best_j);
            dy_recovered(i) = dy_paired(best_j);
            used_new_idx(best_j) = true;
        end
    end
end