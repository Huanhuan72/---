clc ; clear ; close all ;
%% 讀取參數
[ n , ns , qb , ti , tc , tf , Tc , ni , nc , np , nh , nx , wf ] = Parameters() ;  % 載入參數
% ⚠️ 務必確認：Parameters() 裡面的 nx 必須改成 6！

t_f = 0.02:0.02:12;
tsize = size(t_f);
tsize = tsize(2);

%% 優化問題條件 (五次多項式版本)
%======================== 不等式條件( A * x <= B ) =========================
A1_theta = zeros(tsize, 6); A1_v = zeros(tsize, 6); A1_a = zeros(tsize, 6);
A2_theta = zeros(tsize, 6); A2_v = zeros(tsize, 6); A2_a = zeros(tsize, 6);

for i = 1:tsize
    t = t_f(i);
    basis_pos = [1, t, t^2, t^3, t^4, t^5];
    basis_vel = [0, 1, 2*t, 3*t^2, 4*t^3, 5*t^4];
    basis_acc = [0, 0, 2, 6*t, 12*t^2, 20*t^3];
    
    A1_theta(i,:) = basis_pos; A2_theta(i,:) = basis_pos;
    A1_v(i,:) = basis_vel;     A2_v(i,:) = basis_vel;
    A1_a(i,:) = basis_acc;     A2_a(i,:) = basis_acc;
end

A1 = [A1_theta; A1_v; A1_a];
A2 = [A2_theta; A2_v; A2_a];

A_col1 = [A1 ; zeros(1800,6) ; -A1 ; zeros(1800,6)];
A_col2 = [zeros(1800,6) ; A2 ; zeros(1800,6) ; -A2];
A = horzcat(A_col1, A_col2);

the1_max = qb(1,1); the2_max = qb(2,1);
v1_max = qb(1,3);   v2_max = qb(1,3);
a1_max = qb(1,4);   a2_max = qb(1,4);
the1_min = qb(1,2); the2_min = qb(2,2);
B = zeros(7200,1);
B(1:600,1)=the1_max;     B(601:1200,1)=v1_max;    B(1201:1800,1)=a1_max;
B(1801:2400,1)=the2_max; B(2401:3000,1)=v2_max;   B(3001:3600,1)=a2_max;
B(3601:4200,1)=-the1_min;B(4201:4800,1)=v1_max;   B(4801:5400,1)=a1_max;
B(5401:6000,1)=-the2_min;B(6001:6600,1)=v2_max;   B(6601:7200,1)=a2_max;

%======================== 等式條件( Aeq * x = Beq ) ========================
t_start = t_f(1);
t_end = t_f(tsize);

C_pos = [1, t_start, t_start^2, t_start^3, t_start^4, t_start^5] - [1, t_end, t_end^2, t_end^3, t_end^4, t_end^5];
C_vel = [0, 1, 2*t_start, 3*t_start^2, 4*t_start^3, 5*t_start^4] - [0, 1, 2*t_end, 3*t_end^2, 4*t_end^3, 5*t_end^4];
C_acc = [0, 0, 2, 6*t_start, 12*t_start^2, 20*t_start^3] - [0, 0, 2, 6*t_end, 12*t_end^2, 20*t_end^3];

C_base = [C_pos; C_vel; C_acc];
C_col1 = [C_base; zeros(3,6)];
C_col2 = [zeros(3,6); C_base];
Aeq = horzcat(C_col1, C_col2);
Beq = zeros(6,1);

%% 激勵軌跡最佳化
xb = 1 ;  % x0範圍
x0 = 2 * rand( size( A , 2 ) , 1 ) * xb - xb ;  % x 初值
options = optimoptions( 'fmincon' , 'Algorithm' , 'sqp' , 'Display' ,'iter','PlotFcn',{@optimplotfval},'MaxIterations' ,2000 );
[ x , Index ] = fmincon( @(x)OptimizationProblem(x) , x0 , A , B , Aeq , Beq , [] , [] , [] , options ) ;  % 軌跡優化計算
X = reshape( x , nx , n ) ;  % 排列 x

%%
% cmd軌跡圖
the1_f = X(1,1) + X(2,1)*Tc + X(3,1)*Tc.^2 + X(4,1)*Tc.^3 + X(5,1)*Tc.^4 + X(6,1)*Tc.^5;
the2_f = X(1,2) + X(2,2)*Tc + X(3,2)*Tc.^2 + X(4,2)*Tc.^3 + X(5,2)*Tc.^4 + X(6,2)*Tc.^5;

vel1_f = X(2,1) + 2*X(3,1)*Tc + 3*X(4,1)*Tc.^2 + 4*X(5,1)*Tc.^3 + 5*X(6,1)*Tc.^4;
vel2_f = X(2,2) + 2*X(3,2)*Tc + 3*X(4,2)*Tc.^2 + 4*X(5,2)*Tc.^3 + 5*X(6,2)*Tc.^4;

acc1_f = 2*X(3,1) + 6*X(4,1)*Tc + 12*X(5,1)*Tc.^2 + 20*X(6,1)*Tc.^3;
acc2_f = 2*X(3,2) + 6*X(4,2)*Tc + 12*X(5,2)*Tc.^2 + 20*X(6,2)*Tc.^3;

joint1cmd = [the1_f' vel1_f' acc1_f'];
joint2cmd = [the2_f' vel2_f' acc2_f'];

%%
%查看是否超過上下界
figure(2)
subplot(3,2,1); plot(Tc,the1_f); yline(qb(1,1)); yline(qb(1,2)); xlim([0,12]); ylim([-1.5,1.5]); title('Axis1 Position'); xlabel('t(s)'); ylabel('position(rad)');
subplot(3,2,2); plot(Tc,the2_f); yline(qb(2,1)); yline(qb(2,2)); xlim([0,12]); ylim([-2.5,2.5]); title('Axis2 Position'); xlabel('t(s)'); ylabel('position(rad)');
subplot(3,2,3); plot(Tc,vel1_f); xlim([0,12]); ylim([-5.5,5.5]); title('Axis1 Velocity'); xlabel('t(s)'); ylabel('velocity(rad/s)');
subplot(3,2,4); plot(Tc,vel2_f); xlim([0,12]); ylim([-5.5,5.5]); title('Axis2 Velocity'); xlabel('t(s)'); ylabel('velocity(rad/s)');
subplot(3,2,5); plot(Tc,acc1_f); xlim([0,12]); ylim([-10.5,10.5]); title('Axis1 Acceleration'); xlabel('t(s)'); ylabel('acceleration(rad/s^{2})');
subplot(3,2,6); plot(Tc,acc2_f); xlim([0,12]); ylim([-10.5,10.5]); title('Axis2 Acceleration'); xlabel('t(s)'); ylabel('acceleration(rad/s^{2})');

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
fig = figure(3);
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