% 各項參數設定

function [ n , ns , qb , ti , tc , tf , Tc , ni , nc , np , nh , nx , wf ] = Parameters()

%================================= Robot ==================================

n = 2 ;  % 關節數
ns = 8 ;  % 標準參數維度

RatedSpeed = 3000 * 2*pi/60 ;  % 馬達額定轉速 (rad/s) , 單位轉換: (rad/s) = (rpm) * 2*pi/60
GearRatio = [ 50 , 50 ] ;  % 各軸減速比

PosBound = [  90  ,  150  ;
             -90  , -150  ] * 0.9 * ( pi / 180 ) ;  % 機構角度限制 (rad)     
VelBound = ( RatedSpeed ./ GearRatio ) * 0.8 ;  % 機構速度限制 (rad/s)
AccBound = VelBound * 2 ;  % 機構加速度限制: 速度限制*自定倍數 (rad/s^2)         

qb = [ PosBound' , VelBound' , AccBound' ] ;  % 激勵軌跡限制 (角度正、角度負、速度、加速度)

%================================== Time ==================================

ti = 0.02 ;  % 鑑別用取樣時間 (sec)
tc = 0.001 ;  % 控制取樣時間 (sec)
tf = 12 ;  % 基礎週期 = 週期結束時間 (sec)

Tc = tc : tc : tf ;  % 控制時間 (sec)

ni = tf / ti ;  % 單週鑑別用資料數 (samples)
nc = tf / tc ;  % 單週控制資料數 (samples)

%================================ Fourier =================================

np = 10 ;  % 週期數

nh = 5 ;  % 諧波數
nx = 2 * nh + 1 ;  % 單軸目標參數數量
wf = ( 2 * pi ) / tf ;  % 基礎頻率 (rad/s) , 單位轉換: (rad/s) = ( 2 * pi ) * (Hz)

end




