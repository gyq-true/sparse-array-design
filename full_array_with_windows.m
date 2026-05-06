function [U_dB_taylor,U_dB_hamming,U_dB_hanning,U_dB_ChebyShev,V_dB_taylor,V_dB_hamming,V_dB_hanning,V_dB_ChebyShev...
    ] = full_array_with_windows(N_side,k,Space,rho_sll_limit_dB,Axis_scale)
% 满阵加窗（-30db泰勒窗 汉宁窗 汉明窗 切比雪夫窗）比较主瓣宽度、旁瓣电平
%   输入为N_side,k,Space,rho_sll_limit_dB,spacefrequency
%   输出为U_dB_taylor,U_dB_hamming,U_dB_hanning,U_dB_ChebyShev,V_dB_taylor,V_dB_hamming,V_dB_hanning,V_dB_ChebyShev
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
%各种窗
%泰勒窗
win_1d = taylorwin(N_side, 5, rho_sll_limit_dB);    
win_2d = win_1d * win_1d';        
w_full_taylor = win_2d(:);        
%汉明窗
win_ham = hamming(N_side);
w_hamming = (win_ham * win_ham.'); 
w_hamming = w_hamming(:);
%汉宁窗
win_han = hann(N_side); % MATLAB 建议使用 hann 替代 hanning
w_hanning = (win_han * win_han.');
w_hanning = w_hanning(:);
%切比雪夫窗
win_c = chebwin(N_side, abs(rho_sll_limit_dB));
w_cheby = (win_c * win_c.');
w_cheby = w_cheby(:);
%生成方向图
Phase=exp(1j * k * (pos_x * u_flat.' + pos_y * v_flat.'));
AF_taylor=Phase.'*w_full_taylor;
AF_hamming=Phase.'*w_hamming;
AF_hann=Phase.'*w_hanning;
AF_ChebyShev=Phase.'*w_cheby;
%重排
AF_reshape_taylor=reshape(AF_taylor,rows,cols);
AF_reshape_hamming=reshape(AF_hamming,rows,cols);
AF_reshape_hanning=reshape(AF_hann,rows,cols);
AF_reshape_ChebyShev=reshape(AF_ChebyShev,rows,cols);

%u轴方向图切片
[~, idx_v0] = min(abs(v_vec));
U_slice_taylor= AF_reshape_taylor(idx_v0, :);
U_slice_hamming=AF_reshape_hamming(idx_v0, :);
U_slice_hanning=AF_reshape_hanning(idx_v0, :);
U_slice_ChebyShev=AF_reshape_ChebyShev(idx_v0, :);
U_dB_taylor = 20 * log10(abs(U_slice_taylor)/max(abs(U_slice_taylor(:))) ); 
U_dB_hamming = 20 * log10(abs(U_slice_hamming)/max(abs(U_slice_hamming(:))) ); 
U_dB_hanning = 20 * log10(abs(U_slice_hanning)/max(abs(U_slice_hanning(:))) ); 
U_dB_ChebyShev = 20 * log10(abs(U_slice_ChebyShev)/max(abs(U_slice_ChebyShev(:))) ); 
%v轴方向图切片
[~, idx_u0] = min(abs(u_vec));
V_slice_taylor = AF_reshape_taylor(:, idx_u0);
V_slice_hamming = AF_reshape_hamming(:, idx_u0);
V_slice_hanning = AF_reshape_hanning(:, idx_u0);
V_slice_ChebyShev = AF_reshape_ChebyShev(:, idx_u0);
V_dB_taylor = 20 * log10(abs(V_slice_taylor)/max(abs(V_slice_taylor(:))));
V_dB_hamming = 20 * log10(abs(V_slice_hamming)/max(abs(V_slice_hamming(:))));
V_dB_hanning= 20 * log10(abs(V_slice_hanning)/max(abs(V_slice_hanning(:))));
V_dB_ChebyShev= 20 * log10(abs(V_slice_ChebyShev)/max(abs(V_slice_ChebyShev(:))));
end
