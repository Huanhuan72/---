% 優化函數建立(Fourier公式)，輸入為 Fourier 公式正弦與餘弦函數之振幅(目標參數ｘ)，輸出為成本函數(Index)

function [ Index ] = OptimizationProblem( x )

% %============================ 計算優化問題指標 =============================
% 
% %--------------------------------- 條件數 ---------------------------------
% 
% Index = cond( Phi ) ;
% 
% clear Phi ;
[ n , ns , qb , ti , tc , tf , Tc , ni , nc , np , nh , nx , wf ] = Parameters() ;
    samT=0.02;
    phi=zeros(7200,9);
    a=1;
    N=nh;
    d1=0.24;d2=0.24;
    for t=samT:samT:tf
        q1=0;q1d=0;q1dd=0;q2=0;q2d=0;q2dd=0;
        q1=q1+sin(wf*1*t)/wf/1*x(1)-cos(wf*1*t)/wf/1*x(1+1)+x(3);
        q1d=q1d+cos(wf*1*t)*x(1)+sin(wf*1*t)*x(1+1);
        q1dd=q1dd-wf*1*sin(wf*1*t)*x(1)+wf*1*cos(wf*1*t)*x(1+1);
        q2=q2+sin(wf*1*t)/wf/1*x(1+2*N+1)-cos(wf*1*t)/wf/1*x(1+1+2*N+1)+x(14);
        q2d=q2d+cos(wf*1*t)*x(1+2*N+1)+sin(wf*1*t)*x(1+1+2*N+1);
        q2dd=q2dd-wf*1*sin(wf*1*t)*x(1+2*N+1)+wf*1*cos(wf*1*t)*x(1+1+2*N+1);
        for k=2:5
            q1=q1+sin(wf*k*t)/wf/k*x(k+3)-cos(wf*k*t)/wf/k*x(k+4);
            q1d=q1d+cos(wf*k*t)*x(k+3)+sin(wf*k*t)*x(k+4);
            q1dd=q1dd-wf*k*sin(wf*k*t)*x(k+3)+wf*k*cos(wf*k*t)*x(k+4);
            q2=q2+sin(wf*k*t)/wf/k*x(k+14)-cos(wf*k*t)/wf/k*x(k+15);
            q2d=q2d+cos(wf*k*t)*x(k+14)+sin(wf*k*t)*x(k+15);
            q2dd=q2dd-wf*k*sin(wf*k*t)*x(k+14)+wf*k*cos(wf*k*t)*x(k+15);
        end
        phi((a:a+1),:)=[q1dd,q1dd+q2dd,-d1*(sin(q2)*q2d^2+2*q1d*sin(q2)*q1d-2*q1dd*cos(q2)-q2dd*cos(q2)),-d1*(cos(q2)*q2d^2+2*q1d*cos(q2)*q2d+2*q1dd*sin(q2)+q2dd*sin(q2)),0,sign(q1d),0,q1d,0;
                          0,q1dd+q2dd,d1*(sin(q2)*q1d^2+q1dd*cos(q2)),d1*(q1d^2*cos(q2)-q1dd*sin(q2)),q2dd,0,sign(q2d),0,q2d];

        a=a+2;
        
    end
    Index=cond(phi);
    clear phi;
    







