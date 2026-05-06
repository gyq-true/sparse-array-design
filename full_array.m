function [U_dB_full_array,V_dB_full_array] = full_array(N_side,k,Space,Axis_scale)
%满阵不加窗比较一下主瓣宽度和主副瓣比
%输入为：N_side,k,Space,spacefrequency,interp_factor
%输出为：U_dB_full_array,V_dB_full_array
%   和全1的满阵阵元激励进行比较
% =============================================================
%阵元坐标设置
x_vec = (-(N_side/2-1):N_side/2) * Space; % 长度为 32，且包含 0,确保（0,0）处恰好有一个阵元
y_vec = (-(N_side/2-1):N_side/2) * Space;
[X, Y] = meshgrid(x_vec, y_vec);
pos_x = X(:); pos_y = Y(:);
N_total=length(pos_x);     % 更新长度
%空域采样设置
u_vec = linspace(-1, 1, Axis_scale); 
v_vec = linspace(-1, 1, Axis_scale);
[U, V] = meshgrid(u_vec, v_vec);
u_flat = U(:); v_flat = V(:);
[rows, cols] = size(U);
%生成方向图
w=ones(N_total,1);
Phase=exp(1j * k * (pos_x * u_flat.' + pos_y * v_flat.'));
AF_full_array=Phase.'*w;
AF_reshape_full_array=reshape(AF_full_array,rows,cols);
[~, idx_v0] = min(abs(v_vec));
U_Slice_full_array=AF_reshape_full_array(idx_v0, :);
U_dB_full_array=20 * log10(abs(U_Slice_full_array)/max(abs(U_Slice_full_array(:))) );
[~, idx_u0] = min(abs(u_vec));
V_Slice_full_array = AF_reshape_full_array(:, idx_u0);
V_dB_full_array=20 * log10(abs(V_Slice_full_array)/max(abs(V_Slice_full_array(:))) );
end