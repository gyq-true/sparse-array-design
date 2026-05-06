function [hpbw_deg] = calculate_hpbw(u_vec, AF_dB)
% 寻找 HPBW 的辅助逻辑
    % 找到主瓣中心两侧 -3dB 的位置
    idx = find(AF_dB >= -3);
    if isempty(idx), hpbw = nan; 
        fprintf('函数calculate_hpbw执行出错，无法计算3dB波束宽度。')
        return; 
    end
    u_left=u_vec(idx(1));
    u_right=u_vec(idx(end));
    hpbw_deg = rad2deg(asin(max(min(u_right,1),-1)) - asin(max(min(u_left,1),-1)));
end

