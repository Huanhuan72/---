function [c, ceq] = CartesianLimits(x, t_f, wf, nx, n, L1, L2)
    % 將一維的解 x 轉回 nx * n 的矩陣
    X_coeff = reshape(x, nx, n);
    t_check = t_f(1:5:end); 
    tsize = length(t_check);
    
    % 初始化 X, Y 陣列
    X_E = zeros(1, tsize);
    Y_E = zeros(1, tsize);
    
    % 根據傅立葉係數重建時間 t_f 內的角度軌跡
    for i = 1:tsize
        Tc = t_check(i);
        % 這裡帶入你原本計算 the1_f 和 the2_f 的長串公式
        the1 = X_coeff(1,1)*sin(wf*1*Tc)/(wf*1) + X_coeff(4,1)*sin(wf*2*Tc)/(wf*2) + ... % (依此類推，請補齊你原本的 5 階傅立葉展開)
               X_coeff(3,1); 
        the2 = X_coeff(1,2)*sin(wf*1*Tc)/(wf*1) + X_coeff(4,2)*sin(wf*2*Tc)/(wf*2) + ... % (依此類推)
               X_coeff(3,2);
        
        % 順向運動學 (假設 L1 = 0.24, L2 = 0.24)
        X_E(i) = L1 * cos(the1) + L2 * cos(the1 + the2);
        Y_E(i) = L1 * sin(the1) + L2 * sin(the1 + the2);
    end
    
    % --- 設定卡式空間的邊界 (請根據你的機台實際狀況修改) ---
    X_max =  0.4;  % X 軸正向極限
    X_min = -0.4;  % X 軸負向極限
    Y_max =  0.4;  % Y 軸正向極限
    Y_min =  0.1;  % Y 軸負向極限 (例如避免撞到基座或桌面)
    
    % --- 建立不等式條件 (要求 c <= 0) ---
    c_X_max = X_E - X_max;
    c_X_min = X_min - X_E;
    c_Y_max = Y_E - Y_max;
    c_Y_min = Y_min - Y_E;
    
    % 將所有限制條件串接成一個直行向量
    c = [c_X_max, c_X_min, c_Y_max, c_Y_min]';
    
    % 等式條件已在外部 Aeq 處理完畢，這裡留空
    ceq = []; 
end