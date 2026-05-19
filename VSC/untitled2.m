clc;
clear;
close all;

%% 定義開迴路傳遞函數
s = tf('s');
G0 = (s^2 + 3*s + 3) / (s^2*(s+1)*(s+10)*(s+20));

%% 畫根軌跡
figure;
rlocus(G0);
grid on;
title('Root Locus');

%% 由最大超越量求阻尼比
OS = 28;   % 最大超越量 28%
zeta = -log(OS/100) / sqrt(pi^2 + (log(OS/100))^2);

disp(['對應阻尼比 zeta = ', num2str(zeta)]);

%% 加上阻尼比線
sgrid(zeta, []);

%% 取得根軌跡資料
[poles, K] = rlocus(G0);

best_err = inf;
best_K = NaN;
best_pole = NaN;

for i = 1:length(K)
    current_poles = poles(:, i);

    for p = current_poles.'
        % 只考慮左半平面的複數極點
        if imag(p) > 0 && real(p) < 0
            zeta_p = -real(p) / abs(p);
            err = abs(zeta_p - zeta);

            if err < best_err
                best_err = err;
                best_K = K(i);
                best_pole = p;
            end
        end
    end
end

%% 顯示結果
disp(['選到的 K = ', num2str(best_K)]);
disp(['主導極點 = ', num2str(best_pole)]);

%% 建立閉迴路系統
T = feedback(best_K * G0, 1);

%% 步階響應
figure;
step(T);
grid on;
title(['Step Response (K = ', num2str(best_K), ')']);

%% 響應特性
info = stepinfo(T);

disp('Step response information:');
disp(info);
disp(['最大超越量 = ', num2str(info.Overshoot), ' %']);