clc ; clear ; close all ;


%% 讀取參數

[ n , ns , qb , ti , tc , tf , Tc , ni , nc , np , nh , nx , wf ] = Parameters() ;  % 載入參數
% [關節數n, 標準參數維度 ns, 激勵軌跡限制 qb, 鑑別用取樣時間 ti, 基礎週期 tf, 控制時間 t_f, 
% 單週鑑別用資料數 ni, 單週控制資料數 tsize, 週期數 np, 諧波數 nh, 單軸目標參數數量 nx, 基礎頻率 wf]

% 鑑別用的取樣時間和控制用的取樣時間不同 %
t_f = 0.02:0.02:12;
tsize = size(t_f);
tsize = tsize(2);
%% 優化問題條件

%======================== 不等式條件( A * x <= B ) =========================
%                              ( for 空間限制 )

A1_theta = zeros(tsize,11);
A1_v = zeros(tsize,11);
A1_a = zeros(tsize,11);

A2_theta = zeros(tsize,11);
A2_v = zeros(tsize,11);
A2_a = zeros(tsize,11);

for i = 1:tsize
    %theta
    for j = 1:5
        if j == 1
             A1_theta(i,j) = sin(wf*j*t_f(i))/(wf*j);A1_theta(i,j+1) = -cos(wf*j*t_f(i))/(wf*j);A1_theta(i,j+2) = 1;
             A2_theta(i,j) = sin(wf*j*t_f(i))/(wf*j);A2_theta(i,j+1) = -cos(wf*j*t_f(i))/(wf*j);A2_theta(i,j+2) = 1;

            %v
            A1_v(i,j) = cos(wf*j*t_f(i));A1_v(i,j+1)=sin(wf*j*t_f(i));A1_v(i,j+2) = 0;
            A2_v(i,j) = cos(wf*j*t_f(i));A2_v(i,j+1)=sin(wf*j*t_f(i));A2_v(i,j+2) = 0;

            %a
            A1_a(i,j) = -(wf*j)*sin(wf*j*t_f(i));A1_a(i,j+1) = (wf*j)*cos(wf*j*t_f(i));A1_a(i,j+2) = 0;
            A2_a(i,j) = -(wf*j)*sin(wf*j*t_f(i));A2_a(i,j+1) = (wf*j)*cos(wf*j*t_f(i));A2_a(i,j+2) = 0;

        else
            A1_theta(i,2*j) = sin(wf*j*t_f(i))/(wf*j);A1_theta(i,2*j+1) = -cos(wf*j*t_f(i))/(wf*j);    
            A2_theta(i,2*j) = sin(wf*j*t_f(i))/(wf*j);A2_theta(i,2*j+1) = -cos(wf*j*t_f(i))/(wf*j);  
            
            %v
            A1_v(i,2*j) = cos(wf*j*t_f(i));A1_v(i,2*j+1)=sin(wf*j*t_f(i));
            A2_v(i,2*j) = cos(wf*j*t_f(i));A2_v(i,2*j+1)=sin(wf*j*t_f(i));

            %a
            A1_a(i,2*j) = -(wf*j)*sin(wf*j*t_f(i));A1_a(i,2*j+1) = (wf*j)*cos(wf*j*t_f(i));
            A2_a(i,2*j) = -(wf*j)*sin(wf*j*t_f(i));A2_a(i,2*j+1) = (wf*j)*cos(wf*j*t_f(i));
        end
    end
end

A1 = [A1_theta;A1_v;A1_a];
A2 = [A2_theta;A2_v;A2_a];

A_col1 = [A1 ;zeros(1800,11);-A1;zeros(1800,11)];
A_col2 = [zeros(1800,11);A2;zeros(1800,11);-A2];
A = horzcat(A_col1,A_col2);

the1_max = qb(1,1);the2_max = qb(2,1);
v1_max = qb(1,3);v2_max = qb(1,3);
a1_max = qb(1,4);a2_max = qb(1,4);

the1_min = qb(1,2);the2_min = qb(2,2);

B = zeros(7200,1);

B(1:600,1)=the1_max;B(601:1200,1)=v1_max;B(1201:1800,1)=a1_max;
B(1801:2400,1)=the2_max;B(2401:3000,1)=v2_max;B(3001:3600,1)=a2_max;

B(3601:4200,1)=-the1_min;B(4201:4800,1)=v1_max;B(4801:5400,1)=a1_max;
B(5401:6000,1)=-the2_min;B(6001:6600,1)=v2_max;B(6601:7200,1)=a2_max;

%======================== 等式條件( Aeq * x = Beq ) ========================
%                         ( for 末端位置會回到原點 )

C1_theta = zeros(1,11);
C1_v = zeros(2,11);
C1_a = zeros(2,11);

C2_theta = zeros(1,11);
C2_v = zeros(2,11);
C2_a = zeros(2,11);

%初始化1-3個值
%theta
C1_theta(1,1) =  (sin(wf*1*t_f(1))-sin(wf*1*t_f(tsize)))/(wf*1);
C1_theta(1,2) = -(cos(wf*1*t_f(1))+cos(wf*1*t_f(tsize)))/(wf*1);
C1_theta(1,3) = 0;
%C1,v
C1_v(1,1) = cos(wf*1*t_f(1));C1_v(1,2) = sin(wf*1*t_f(1));C1_v(1,3) = 0;
C1_v(2,1) = cos(wf*1*t_f(tsize));C1_v(2,2) = sin(wf*1*t_f(tsize));C1_v(2,3) = 0;
%C1,a
C1_a(1,1) = -(wf*1)*sin(wf*1*t_f(1));C1_a(1,2)=(wf*1)*cos(wf*1*t_f(1));C1_a(1,3) = 0;
C1_a(2,1) = -(wf*1)*sin(wf*1*t_f(tsize));C1_a(2,2)=(wf*1)*cos(wf*1*t_f(tsize));C1_a(2,3) = 0;
%C2,v
C2_v(1,1) = cos(wf*1*t_f(1));C2_v(1,2) = sin(wf*1*t_f(1));C2_v(1,3) = 0;
C2_v(2,1) = cos(wf*1*t_f(tsize));C2_v(2,2) = sin(wf*1*t_f(tsize));C2_v(2,3) = 0;
%C2,a
C2_a(1,1) = -(wf*1)*sin(wf*1*t_f(1));C2_a(1,2)=(wf*1)*cos(wf*1*t_f(1));C2_a(1,3) = 0;
C2_a(2,1) = -(wf*1)*sin(wf*1*t_f(tsize));C2_a(2,2)=(wf*1)*cos(wf*1*t_f(tsize));C2_a(2,3) = 0;

for i = 2:5
    %theta
    C1_theta(1,2*i) =  (sin(wf*i*t_f(1))-sin(wf*i*t_f(tsize)))/(wf*i);
    C1_theta(1,2*i+1) = -(cos(wf*i*t_f(1))+cos(wf*i*t_f(tsize)))/(wf*i);
    
    C2_theta(1,2*i) =  (sin(wf*i*t_f(1))-sin(wf*i*t_f(tsize)))/(wf*i);
    C2_theta(1,2*i+1) = -(cos(wf*i*t_f(1))+cos(wf*i*t_f(tsize)))/(wf*i);

    %C1 v,a
    C1_v(1,2*i) = cos(wf*i*t_f(1));C1_v(1,2*i+1) = sin(wf*i*t_f(1));
    C1_v(2,2*i) = cos(wf*i*t_f(tsize));C1_v(2,2*i+1) = sin(wf*i*t_f(tsize));

    C1_a(1,2*i) = -(wf*i)*sin(wf*i*t_f(1));C1_a(1,2*i)=(wf*i)*cos(wf*i*t_f(1));
    C1_a(2,2*i) = -(wf*i)*sin(wf*i*t_f(tsize));C1_a(2,2*i)=(wf*i)*cos(wf*i*t_f(tsize));
    
    %C2 v,a
    C2_v(1,2*i) = cos(wf*i*t_f(1));C2_v(1,2*i+1) = sin(wf*i*t_f(1));
    C2_v(2,2*i) = cos(wf*i*t_f(tsize));C2_v(2,2*i+1) = sin(wf*i*t_f(tsize));

    C1_a(1,2*i) = -(wf*i)*sin(wf*i*t_f(1));C1_a(1,2*i+1)=(wf*2)*cos(wf*i*t_f(1));
    C2_a(2,2*i) = -(wf*i)*sin(wf*i*t_f(tsize));C2_a(2,2*i+1)=(wf*2)*cos(wf*i*t_f(tsize));

end

C1 = [C1_theta;C1_v;C1_a];
C2 = [C2_theta;C2_v;C2_a];

C_col1 = [C1;zeros(5,11)];
C_col2 = [zeros(5,11);C2];
Aeq = horzcat(C_col1,C_col2);

Beq = zeros(10,1);


%% 激勵軌跡最佳化 (加入卡式空間事後驗證機制)
fprintf('準備開始 PSO 最佳化...\n');

% 1. 設定變數數量與搜尋上下界
nvars = size(A, 2); % 傅立葉係數的總數量
xb = 1;             % 你原本設定的係數搜尋範圍
lb = -xb * ones(nvars, 1); % 下界
ub =  xb * ones(nvars, 1); % 上界
% 2. 設定 PSO 選項
% SwarmSize 決定粒子數量，越大找得越準但越慢
options = optimoptions('particleswarm', ...
    'SwarmSize', 150, ...
    'MaxIterations', 100, ...
    'Display', 'iter');
% 3. 執行 PSO
% 注意：這裡把 A, B, Aeq, Beq, nx, n, wf, Tc 都傳進去給 Fitness Function 評分
[x_pso, fval, exitflag, output] = particleswarm(...
    @(x) PSO_FitnessFunction(x, A, B, Aeq, Beq, nx, n, wf, Tc), ...
    nvars, lb, ub, options);

fprintf('>>> PSO 最佳化完成！最終適應度分數: %.4f\n\n', fval);

% 4. 取得最佳解並還原成 X 矩陣
x = x_pso(:); % 確保是行向量
X = reshape(x, nx, n);

% ==========================================================
% 下方的程式碼維持不變 (計算 the1_f, the2_f, vel1_f 等等...)

the1_f =  X(1,1)*sin(wf *1 *Tc)/(wf *1)+X(4,1)*sin(wf *2 *Tc)/(wf *2) + X(6,1)*sin(wf *3 *Tc)/(wf *3) + ...
               X(8,1)*sin(wf *4 *Tc)/(wf *4) + X(10,1)*sin(wf *5 *Tc)/(wf *5) + X(2,1)*(-cos(wf *1 *Tc)/(wf *1))+ ...
               X(5,1)*(-cos(wf *2 *Tc)/(wf *2)) + X(7,1)*(-cos(wf *3 *Tc)/(wf *3))+ X(9,1)*(-cos(wf *4 *Tc)/(wf *4))+ ...
               X(11,1)*(-cos(wf *5 *Tc)/(wf *5)) + X(3,1);
    the2_f =  X(1,2)*sin(wf *1 *Tc)/(wf *1)+X(4,2)*sin(wf *2 *Tc)/(wf *2) + X(6,2)*sin(wf *3 *Tc)/(wf *3) + ...
               X(8,2)*sin(wf *4 *Tc)/(wf *4) + X(10,2)*sin(wf *5 *Tc)/(wf *5) + X(2,2)*(-cos(wf *1 *Tc)/(wf *1))+ ...
               X(5,2)*(-cos(wf *2 *Tc)/(wf *2)) + X(7,2)*(-cos(wf *3 *Tc)/(wf *3))+ X(9,2)*(-cos(wf *4 *Tc)/(wf *4))+ ...
               X(11,2)*(-cos(wf *5 *Tc)/(wf *5)) + X(3,2);
vel1_f = X(1,1)* cos((wf *1 ) * Tc) + X(4,1)* cos((wf *2 )* Tc) + X(6,1)* cos((wf* 3) * Tc) + X(8,1)* cos((wf* 4) * Tc)+ ...
           X(10,1)* cos((wf* 5) * Tc) + X(2,1)* sin((wf *1 )* Tc) + X(5,1)* sin((wf *2 )* Tc) + X(7,1)* sin((wf *3 )* Tc)+ ...
           X(9,1)* sin((wf *4 )* Tc) + X(11,1)* sin((wf *5 )* Tc);

vel2_f = X(1,2)* cos((wf *1 ) * Tc) + X(4,2)* cos((wf *2 )* Tc) + X(6,2)* cos((wf* 3) * Tc) + X(8,2)* cos((wf* 4) * Tc)+ ...
           X(10,2)* cos((wf* 5) * Tc) + X(2,2)* sin((wf *1 )* Tc) + X(5,2)* sin((wf *2 )* Tc) + X(7,2)* sin((wf *3 )* Tc)+ ...
           X(9,2)* sin((wf *4 )* Tc) + X(11,2)* sin((wf *5 )* Tc);

acc1_f = X(1,1) *wf *1 * -sin((wf *1 ) * Tc) + X(4,1)*wf *2 * -sin((wf *2 )* Tc) + X(6,1)*wf *3 * -sin((wf* 3) * Tc) + X(8,1)*wf *4 * -sin((wf* 4) * Tc)+ ...
           X(10,1)*wf *5 * -sin((wf* 5) * Tc) + X(2,1)*wf *1 *cos((wf *1 )* Tc) + X(5,1)*wf *2* cos((wf *2 )* Tc) + X(7,1)*wf *3* cos((wf *3 )* Tc)+ ...
           X(9,1)*wf *4* cos((wf *4 )* Tc) + X(11,1)*wf *5* cos((wf *5 )* Tc);

acc2_f = X(1,2) *wf *1 * -sin((wf *1 ) * Tc) + X(4,2)*wf *2 * -sin((wf *2 )* Tc) + X(6,2)*wf *3 * -sin((wf* 3) * Tc) + X(8,2)*wf *4 * -sin((wf* 4) * Tc)+ ...
           X(10,2)*wf *5 * -sin((wf* 5) * Tc) + X(2,2)*wf *1 *cos((wf *1 )* Tc) + X(5,2)*wf *2* cos((wf *2 )* Tc) + X(7,2)*wf *3* cos((wf *3 )* Tc)+ ...
           X(9,2)*wf *4* cos((wf *4 )* Tc) + X(11,2)*wf *5* cos((wf *5 )* Tc);
% 
joint1cmd = [the1_f' vel1_f' acc1_f'];
joint2cmd = [the2_f' vel2_f' acc2_f'];

%%
%查看是否超過上下界
figure(2)
subplot(3,2,1)
plot(Tc,the1_f);
yline(qb(1,1));
yline(qb(1,2));
xlim([0,12]);
ylim([-1.5,1.5]);
title('Axis1 Position');
xlabel('t(s)');
ylabel('position(rad)');

subplot(3,2,2)
plot(Tc,the2_f);
yline(qb(2,1));
yline(qb(2,2));
xlim([0,12]);
ylim([-2.5,2.5]);
title('Axis2 Position');
xlabel('t(s)');
ylabel('position(rad)');


subplot(3,2,3)
plot(Tc,vel1_f);
xlim([0,12]);
ylim([-5.5,5.5]);
title('Axis1 Velocity');
xlabel('t(s)');
ylabel('velocity(rad/s)');

subplot(3,2,4)
plot(Tc,vel2_f);
xlim([0,12])
ylim([-5.5,5.5])
title('Axis2 Velocity');
xlabel('t(s)');
ylabel('velocity(rad/s)');

subplot(3,2,5)
plot(Tc,acc1_f);
xlim([0,12])
ylim([-10.5,10.5]);
title('Axis1 Acceleration');
xlabel('t(s)');
ylabel('acceleration(rad/s^{2})');

subplot(3,2,6)
plot(Tc,acc2_f);
xlim([0,12])
ylim([-10.5,10.5])
title('Axis2 Acceleration');
xlabel('t(s)');
ylabel('acceleration(rad/s^{2})');

%%
%順向運動學
tf_s = size(Tc);
T_01 = zeros(4);
T_12 = zeros(4);
X1_E = zeros(2,tf_s(2));
X2_E = zeros(2,tf_s(2));

for i_f = 1:tf_s(2)
    T_01 = [cos(the1_f(i_f)) -cos(0)*sin(the1_f(i_f))  sin(0)*sin(the1_f(i_f)) 0.24*cos(the1_f(i_f));
            sin(the1_f(i_f))  cos(0)*cos(the1_f(i_f)) -sin(0)*cos(the1_f(i_f)) 0.24*sin(the1_f(i_f));
            0             sin(0)               cos(0)              0;
            0             0                    0                   1];

    T_12 = [cos(the2_f(i_f)) -cos(0)*sin(the2_f(i_f))  sin(0)*sin(the2_f(i_f)) 0.24*cos(the2_f(i_f));
            sin(the2_f(i_f))  cos(0)*cos(the2_f(i_f)) -sin(0)*cos(the2_f(i_f)) 0.24*sin(the2_f(i_f));
            0             sin(0)               cos(0)              0;
            0             0                    0                   1];

    X1_E_moment = T_01*transpose([0 0 0 1]);
    X2_E_moment = T_01*T_12*transpose([0 0 0 1]);
    
    X1_E(1,i_f) = X1_E_moment(1,1);
    X1_E(2,i_f) = X1_E_moment(2,1);
    X2_E(1,i_f) = X2_E_moment(1,1);
    X2_E(2,i_f) = X2_E_moment(2,1);
end

the1_f = repmat(the1_f',10,1);the2_f = repmat(the2_f',10,1);
vel1_f = repmat(vel1_f',10,1);vel2_f = repmat(vel2_f',10,1);
acc1_f = repmat(acc1_f',10,1);acc2_f = repmat(acc2_f',10,1);

%%
%模擬動畫
fig = figure(3)
h = animatedline;
axis([-0.5,0.5,-0.5,0.5]);
base = [0 0];
pic_num = 1;

xlabel('X','Fontsize',10);
ylabel('Y','Fontsize',10);
for k = 1:100:tf_s(2)
    caption = sprintf('time %d', Tc(k));
    title(caption,'Fontsize',15);
    link1 = line([0,X1_E(1,k)],[0,X1_E(2,k)],'Marker','o');
    link2 = line([X1_E(1,k),X2_E(1,k)],[X1_E(2,k),X2_E(2,k)],'Marker','o');
    link = [link1,link2];
    addpoints(h,X2_E(1,k),X2_E(2,k));
    drawnow
    F = getframe(fig);
    I = frame2im(F);
    [I,map] = rgb2ind(I,256);
    if pic_num ==1
        imwrite(I,map,'path.gif','gif','Loopcount',inf,'DelayTime',0.05);
    else
        imwrite(I,map,'path.gif','gif','WriteMode','append','DelayTime',0.05);
    end
    pic_num = pic_num+1;
    pause(0.05)
    delete(link)

end


%%
%儲存txt
fid1 = fopen('joint.txt','wt');
[m,n] = size(joint1cmd);
for i = 1:1:m*10
    fprintf(fid1,'%f\t%f\t%f\t%f\t%f\t%f\t\n',the1_f(i),the2_f(i),vel1_f(i),vel2_f(i),acc1_f(i),acc2_f(i));
end
fclose(fid1);
