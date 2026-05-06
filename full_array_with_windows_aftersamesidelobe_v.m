function [Slice_V_taylor,Slice_V_ChebyShev] =  full_array_with_windows_aftersamesidelobe_v(N_side,k,Space,max_v_sll_val,Axis_scale)
% v轴满阵加窗（泰勒窗 切比雪夫窗）在相同旁瓣电平下比较主瓣宽度
%   此处提供详细说明
% =============================================================
%阵元坐标设置
x_vec = (-(N_side/2-1):N_side/2) * Space; % 长度为 32，且包含 0,确保（0,0）处恰好有一个阵元
y_vec = (-(N_side/2-1):N_side/2) * Space;
[X, Y] = meshgrid(x_vec, y_vec);
pos_x = X(:); pos_y = Y(:);
N_total=length(pos_x);     % 更新长度
%空域采样设置
u_vec = linspace(-1, 1, Axis_scale); %角域1度的分辨率等效到u/v域，则u/v域约为115个网格点
v_vec = linspace(-1, 1, Axis_scale);
[U, V] = meshgrid(u_vec, v_vec);
u_flat = U(:); v_flat = V(:);
[rows, cols] = size(U);
%各种窗
%泰勒窗
win_1d = taylorwin(N_side, 5, max_v_sll_val);    
win_2d = win_1d * win_1d';        
w_full_taylor = win_2d(:);        
%切比雪夫窗
win_c = chebwin(N_side, abs(max_v_sll_val));
w_cheby = (win_c * win_c.');
w_cheby = w_cheby(:);
%生成方向图
Phase=exp(1j * k * (pos_x * u_flat.' + pos_y * v_flat.'));
AF_taylor=Phase.'*w_full_taylor;
AF_ChebyShev=Phase.'*w_cheby;
%重排
AF_reshape_taylor=reshape(AF_taylor,rows,cols);
AF_reshape_ChebyShev=reshape(AF_ChebyShev,rows,cols);
%v轴方向图切片
[~, idx_u0] = min(abs(u_vec));
V_slice_taylor = AF_reshape_taylor(:, idx_u0);
V_slice_ChebyShev = AF_reshape_ChebyShev(:, idx_u0);
Slice_V_taylor = 20 * log10(abs(V_slice_taylor)/max(abs(V_slice_taylor(:))));
Slice_V_ChebyShev= 20 * log10(abs(V_slice_ChebyShev)/max(abs(V_slice_ChebyShev(:))));
end
