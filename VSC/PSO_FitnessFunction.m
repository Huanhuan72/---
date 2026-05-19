function cost = PSO_FitnessFunction(x, A, B, Aeq, Beq, nx, n, wf, Tc)
    % 確保 x 是行向量 (PSO 傳進來通常是列向量 1xN)
    x = x(:); 
    
    % ---------------------------------------------------
    % 1. 目標函數 (Base Cost)
    % 這裡呼叫你原本的 OptimizationProblem，或者是計算條件數
    % base_cost = OptimizationProblem(x); 
    % 為了示範，若你沒有提供原始的 OptimizationProblem，這裡先設為 0
    % 實際上你要把希望「最小化」的目標值放在這裡
    base_cost = 0; 
    
    % ---------------------------------------------------
    % 2. 懲罰權重設定
    % 權重必須極大，才能有效嚇阻粒子進入違規區
    penalty = 0;
    W_ineq = 1e6;  % 超出關節速度/加速度極限的懲罰
    W_eq   = 1e6;  % 不滿足週期性的懲罰
    W_cart = 1e6;  % 撞到卡氏空間邊界的懲罰
    
    % ---------------------------------------------------
    % 3. 檢查不等式條件 (A*x <= B)
    ineq_check = A * x - B;
    if any(ineq_check > 0)
        % 計算超出量的總和，乘上權重
        penalty = penalty + sum(ineq_check(ineq_check > 0)) * W_ineq;
    end
    
    % 4. 檢查等式條件 (Aeq*x == Beq)
    eq_check = abs(Aeq * x - Beq);
    % 容許 1e-4 的數值計算誤差
    if sum(eq_check) > 1e-4 
        penalty = penalty + sum(eq_check) * W_eq;
    end
    
    % ---------------------------------------------------
    % 5. 檢查卡氏空間限制 (防撞)
    % 先將 x 還原成兩軸位置 (使用你原本的算式)
    X_mat = reshape(x, nx, n);
    
    the1_f = X_mat(1,1)*sin(wf*1*Tc)/(wf*1) + X_mat(4,1)*sin(wf*2*Tc)/(wf*2) + X_mat(6,1)*sin(wf*3*Tc)/(wf*3) + ...
             X_mat(8,1)*sin(wf*4*Tc)/(wf*4) + X_mat(10,1)*sin(wf*5*Tc)/(wf*5) + X_mat(2,1)*(-cos(wf*1*Tc)/(wf*1))+ ...
             X_mat(5,1)*(-cos(wf*2*Tc)/(wf*2)) + X_mat(7,1)*(-cos(wf*3*Tc)/(wf*3))+ X_mat(9,1)*(-cos(wf*4*Tc)/(wf*4))+ ...
             X_mat(11,1)*(-cos(wf*5*Tc)/(wf*5)) + X_mat(3,1);
             
    the2_f = X_mat(1,2)*sin(wf*1*Tc)/(wf*1) + X_mat(4,2)*sin(wf*2*Tc)/(wf*2) + X_mat(6,2)*sin(wf*3*Tc)/(wf*3) + ...
             X_mat(8,2)*sin(wf*4*Tc)/(wf*4) + X_mat(10,2)*sin(wf*5*Tc)/(wf*5) + X_mat(2,2)*(-cos(wf*1*Tc)/(wf*1))+ ...
             X_mat(5,2)*(-cos(wf*2*Tc)/(wf*2)) + X_mat(7,2)*(-cos(wf*3*Tc)/(wf*3))+ X_mat(9,2)*(-cos(wf*4*Tc)/(wf*4))+ ...
             X_mat(11,2)*(-cos(wf*5*Tc)/(wf*5)) + X_mat(3,2);
             
    % 順向運動學
    L1 = 0.24; L2 = 0.24;
    X_E = L1 * cos(the1_f) + L2 * cos(the1_f + the2_f);
    Y_E = L1 * sin(the1_f) + L2 * sin(the1_f + the2_f);
    X_elbow = L1 * cos(the1_f);
    Y_elbow = L1 * sin(the1_f);
    
    % 設定安全邊界 (與你原本設定相同)
    X_max_limit = 1;
    X_min_limit = 0;
    Y_max_limit = 0.4;
    Y_min_limit = -0.4;
    
    % 檢查末端與手肘是否有超出邊界
    all_X = [X_E, X_elbow];
    all_Y = [Y_E, Y_elbow];
    
    out_X_max = sum(max(0, all_X - X_max_limit));
    out_X_min = sum(max(0, X_min_limit - all_X));
    out_Y_max = sum(max(0, all_Y - Y_max_limit));
    out_Y_min = sum(max(0, Y_min_limit - all_Y));
    
    total_cart_violation = out_X_max + out_X_min + out_Y_max + out_Y_min;
    
    if total_cart_violation > 0
        penalty = penalty + total_cart_violation * W_cart;
    end
    
    % ---------------------------------------------------
    % 6. 輸出最終分數給 PSO 評估
    cost = base_cost + penalty;
end