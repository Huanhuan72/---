%% 實驗練習 1_4：步階響應與內部訊號分析 (M-file)

% 1. 定義時間向量 (30秒)
t = 0:0.01:30;

% 2. 建立系統轉移函數
% G1(s) = 1 / (s^2 - 2s + 1)
G1 = tf(1, [1 -2 1]);

% G2(s) = 12(s-1)^2 / ((s+3)(s+1)^2(s+2))
% 使用 zpk (零點、極點、增益) 來定義 G2 比較直觀且不會出錯
z = [1; 1];
p = [-3; -1; -1; -2];
k = 12;
G2 = zpk(z, p, k);

% 3. 求取系統輸出 y 的閉迴路系統 (R 到 Y)
% 理論上 T_y = G1*G2 / (1 + G1*G2)
% 使用 minreal 來消除數值計算上的不穩定對消項，觀察純外部表現
G_open = minreal(G1 * G2);
sys_y = feedback(G_open, 1);

% 4. 求取內部訊號 u 的閉迴路系統 (R 到 U)
% T_u = G1 / (1 + G1*G2)
% 這裡不使用 minreal，還原真實世界中無法完美對消導致的內部發散狀況
sys_u = feedback(G1, G2);

% 5. 繪製結果
figure('Name', '實驗練習 1_4：步階響應');

% 繪製 y 的響應 (表現為穩定)
subplot(2,1,1);
step(sys_y, t);
title('系統輸出 y 的步階響應 (外部穩定)');
xlabel('Time (seconds)'); ylabel('Amplitude (y)');
grid on;

% 繪製 u 的響應 (內部發散)
subplot(2,1,2);
step(sys_u, t);
title('內部訊號 u 的步階響應 (內部不穩定，訊號發散)');
xlabel('Time (seconds)'); ylabel('Amplitude (u)');
grid on;