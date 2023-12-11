clc;
clear;
close all;

tic;            %Start timer(計時耗費時間用，Command Window跳出Elapsed time才是程式結束)
                %程式跑的時間約為40秒，根據每台電腦的性能或是angle解析度有所差異(d_ang越小時間越長，反之則時間越短)
%% Initialization
d_ang = 0.01;                           %angle解析度 (-pi到pi，每隔多少取一次值)
fs = 10;                                %取樣率(DTFT使用)
ts = 1/fs;                              %取樣週期(DTFT使用)
ang = -pi-d_ang:d_ang:pi+d_ang;         %angle陣列
length_of_cos = 50;                     %每種cos的取樣個數(homework範例為50)
cos_arr1=zeros(1,length_of_cos);        %0.2pi的cos
cos_arr2=zeros(1,length_of_cos);        %0.4pi的cos
cos_arr3=zeros(1,length_of_cos);        %0.8pi的cos
thetas_c = zeros(1,10);                 
thetas_x = zeros(1,10);
aa = 0.98;                              %H1分子的系數，測試一些東西可以直接改(預設為0.98)
%% Generate message start

for i=1:length_of_cos
    cos_arr1(i)=cos(0.2*pi*i);
    cos_arr2(i)=cos(0.4*pi*i);
    cos_arr3(i)=cos(0.8*pi*i);
end

ham_window=zeros(1,length_of_cos+1);

for i=1:length_of_cos+1                                   %Hamming
    ham_window(i) = 0.54-0.46*cos(2*pi*(i/length_of_cos-1));
end

fin = zeros(1,length_of_cos);
fin2 = zeros(1,length_of_cos);
fin3 = zeros(1,length_of_cos);


for i=1:length_of_cos                                   %Hamming
    fin(i)=cos_arr1(i)*ham_window(i);
    fin2(i)=cos_arr2(i)*ham_window(i);
    fin3(i)=cos_arr3(i)*ham_window(i);
end

z_p = zeros(1,length(ang)-length_of_cos*3);             %zero-padding
c = cat(2,fin3,fin,fin2,z_p);                           %3 cosine assembling


%% Plot message start

[X1,x1,df1] = fftseq(c,ts,d_ang);                       %DTFT
f = ((0:df1:df1*(length(x1)-1))-fs/2)*0.2*pi;

figure(1);                                              %Plot figure c in homwork
subplot(2,1,1);
plot(c);
title('Waveform of signal x[n]');
grid on
xlabel("Sample number (n)")
subplot(2,1,2);
plot(f,fftshift(abs(X1)));
title('Magnitude of DTFT of x[n]');
grid on
xlabel("\omega");
xticks([-pi -0.8*pi -0.6*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.6*pi 0.8*pi pi])
xticklabels({'-\pi','-0.8\pi','-0.6\pi','-0.4\pi','-0.2\pi','0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
ylabel("|X(e^{j\omega})|");

%% Chanel 1 (H1) start
theta1 = exp(1i*0.8*pi);                %分子(分數上半)角度                    
theta2 = exp(1i*0.4*pi);                %分母(分數下半)角度

thetas_c(1) = aa*theta1;                %紀錄角度，以畫出pole-zero unit circle的圖形             
thetas_c(2) = conj(thetas_c(1));        %尾綴為c的是畫出zero，尾綴為x則是畫出pole
thetas_x(1) = 0.8*theta2;
thetas_x(2) = conj(thetas_x(1));

H1u_arr = zeros(1,length(ang));         %紀錄分子(分數上半，故尾綴為u)
H1d_arr = zeros(1,length(ang));         %紀錄分母(分數下半，故尾綴為d)
H1 = zeros(1,length(ang));              %紀錄分子除分母的最終值
pH1 = zeros(1,length(ang));             %紀錄分子除分母的最終值的相位

for k = 1:length(ang)                   %帶入訊號x[n]並計算出經過第一個channel(H1)後的狀態
    z = power(exp(1i*ang(k)),-1);
    H1u_arr(k)= (1-aa*theta1*z)*(1-aa*conj(theta1)*z);
    H1d_arr(k) = (1-0.8*theta2*z)*(1-0.8*conj(theta2)*z);
    H1(k) = H1u_arr(k)/H1d_arr(k);
    pH1(k) = atan(imag(H1(k))/real(H1(k)));
end

%% Chanel 2 (H2) start
theta3 = exp(1i*0.15*pi);               %分子分母角度一樣

H2u_arr = zeros(1,length(ang));         %紀錄分子(分數上半，故尾綴為u)
H2d_arr = zeros(1,length(ang));         %紀錄分母(分數下半，故尾綴為d)
H2_t = zeros(1,length(ang));            %由於第二個channel有平方，經過第一次方須先暫存
H2 = zeros(1,length(ang));              %紀錄H2最終值
pH2 = zeros(1,length(ang));             %紀錄H2最終值的相位

for m = 1:4                                 %帶入訊號x[n]*H1並計算出經過第二個channel(H1)後的狀態
                                            %m = 1,2,3,4
    theta3 = exp(1i*(0.15*pi+0.02*pi*m));   %角度隨著m改變

    thetas_c(2*m+1) = (1/0.95)*theta3;      %紀錄角度，以畫出pole-zero unit circle的圖形
    thetas_c(2*m+2) = (1/0.95)*conj(theta3);%尾綴為c的是畫出zero，尾綴為x則是畫出pole
    thetas_x(2*m+1) = 0.95*theta3;
    thetas_x(2*m+2) = 0.95*conj(theta3);

    for n = 1:length(ang)                   %經過第一次方
        z = power(exp(1i*ang(n)),-1);
        H2u_arr(n) = (0.95*conj(theta3)-z)*(0.95*theta3-z);
        H2d_arr(n) = (1-0.95*theta3*z)*(1-0.95*conj(theta3)*z);
        H2_t(n) = power((H2u_arr(n)/H2d_arr(n)),2);         %經過第二次方(這邊直接用power(*,2)運算)
    end
    
    if (m==1)                               %將H2作連乘
        H2 = H2_t;
    else
        H2 = H2.*H2_t;
    end
end

H = H1.*H2;                                 %經過兩個通道
pH = zeros(1,length(H));                    %紀錄兩個通道的相位

%% Unwrapped phase start

for l = 1:length(H)                         %算出相位(arctan(Im(*)/Re(*)))
    pH(l) = atan2(imag(H(l)),real(H(l)));
end
lmax = islocalmax(pH);                      %找出區域最大值
lmax = find(lmax);

figure(2);                                  %Plot figure a(1) in homework
subplot(2,1,1);
plot(ang,pH);
title('Principle Value of Phase Response');
grid on
xlabel("\omega");
xticks([-pi -0.8*pi -0.6*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.6*pi 0.8*pi pi])
xticklabels({'-\pi','-0.8\pi','-0.6\pi','-0.4\pi','-0.2\pi','0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
ylabel("ARG[H(e^{j\omega})]");


for k = 2:(length(lmax)-1)                  %Unwrap
    for l = lmax(k):length(pH)              %只要遇到區域最大值，便將後面的項數值全部減2pi      
        pH(l) = pH(l)-2*pi;
    end
end

diff = abs(pH(round(length(pH)/2)));        %查看整個pH的中間項數數值

for l = 1:length(pH)                        %調整至0位點
    pH(l)= pH(l)+diff;
end


subplot(2,1,2);                             %Plot figure a(2) in homework
plot(ang,pH);
title('Unwrapped Phase Response');
grid on
xlabel("\omega");
xticks([-pi -0.8*pi -0.6*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.6*pi 0.8*pi pi])
xticklabels({'-\pi','-0.8\pi','-0.6\pi','-0.4\pi','-0.2\pi','0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
ylabel("arg[H(e^{j\omega})]");

%% Group delay start
slope = zeros(1,length(H));                 %Group delay為對phase做微分並冠負號
for k = 1:length(pH)-1                      %取斜率(此刻項數值-下一刻項數值)/angle解析度
    slope(k) = -(pH(k+1)-pH(k))/d_ang;
end


figure(3);                                  %Plot figure b in homework
subplot(2,1,1);
plot(ang,slope);
title('Group delay of H(z)');
grid on
xlabel("\omega");
xticks([-pi -0.8*pi -0.6*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.6*pi 0.8*pi pi])
xticklabels({'-\pi','-0.8\pi','-0.6\pi','-0.4\pi','-0.2\pi','0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
ylabel("grd[H(e^{j\omega})]");
subplot(2,1,2);
plot(ang,abs(H));
title('Magnitude of Frequency Response');
grid on
xlabel("\omega");
xticks([-pi -0.8*pi -0.6*pi -0.4*pi -0.2*pi 0 0.2*pi 0.4*pi 0.6*pi 0.8*pi pi])
xticklabels({'-\pi','-0.8\pi','-0.6\pi','-0.4\pi','-0.2\pi','0','0.2\pi','0.4\pi','0.6\pi','0.8\pi','\pi'})
ylabel("[H(e^{j\omega})]");

%% Channel 1 LCCDE start
syms z                                          %以z做為多項式
H1u= (1-aa*theta1*z)*(1-aa*conj(theta1)*z);     %分子(後綴為u)
H1d = (1-0.8*theta2*z)*(1-0.8*conj(theta2)*z);  %分母(後綴為d)
cx1 = coeffs(H1u);                              %取出分子係數(做為LCCDE x 的系數)
cy1 = coeffs(H1d);                              %取出分母係數(做為LCCDE y 的系數)

c_h1_temp = zeros(1,length(c));                 %暫存經過channel 1的結果值

for i=1:length(c)                               %根據LCCDE做計算  
    if(i==1)
        c_h1_temp(i) = c(i);
    elseif(i==2)
        c_h1_temp(i) = -1*(cy1(2))*c_h1_temp(i-1)+c(i)+(cx1(2))*c(i-1);
    else
        c_h1_temp(i) = -1*(cy1(2))*c_h1_temp(i-1)-1*(cy1(3))*c_h1_temp(i-2)+c(i)+(cx1(2))*c(i-1)+(cx1(3))*c(i-2);
    end
end

figure(4);                                      %這張圖是多畫的，目的是為了要看原本訊號經過channel 1是否有效果
subplot(2,1,1);
plot(c);
title('Waveform of signal x[n]');
grid on
xlabel("Sample number (n)");
subplot(2,1,2);
plot(c_h1_temp);
title('Waveform of signal x[n] after H1');
grid on
xlabel("Sample number (n)");

%% Channel 2 LCCDE start
syms z                                              %以z做為多項式
c_h2_temp = zeros(1,length(c));                     %暫存
c_h2_t = zeros(1,length(c));                        %暫存
c_h2 = zeros(1,length(c));                          %儲存結果值


for i = 1:length(c_h1_temp)
    c_h2_t(i) = c_h1_temp(i);                       %將channel 1的結果值移進暫存
end

for m = 1:4                                         %m=1,2,3,4
    theta3 = exp(1i*(0.15*pi+0.02*pi*m));           %角度隨著m改變    
    H2u= (0.95*conj(theta3)-z)*(0.95*theta3-z);     %分子(後綴為u)
    H2d= (1-0.95*theta3*z)*(1-0.95*conj(theta3)*z); %分母(後綴為d)
    cx2 = coeffs(H2u);                              %取出分子係數(做為LCCDE x 的系數)
    cy2 = coeffs(H2d);                              %取出分母係數(做為LCCDE y 的系數)    
    for i=1:length(c)                               %根據LCCDE做計算
        if(i==1)
            c_h2_temp(i) = (cx2(1))*c_h2_t(i);
        elseif(i==2)
            c_h2_temp(i) = -1*(cy2(2))*c_h2_temp(i-1)+(cx2(1))*c_h2_t(i)+(cx2(2))*c_h2_t(i-1);
        else
            c_h2_temp(i) = -1*(cy2(2))*c_h2_temp(i-1)-1*(cy2(3))*c_h2_temp(i-2)+(cx2(1))*c_h2_t(i)+(cx2(2))*c_h2_t(i-1)+(cx2(3))*c_h2_t(i-2);
        end
    end
    for i=1:length(c)                               %由於有平方(二次方)，所以再次計算
        if(i==1)
            c_h2(i) = (cx2(1))*c_h2_temp(i);
        elseif(i==2)
            c_h2(i) = -1*(cy2(2))*c_h2(i-1)+(cx2(1))*c_h2_temp(i)+(cx2(2))*c_h2_temp(i-1);
        else
            c_h2(i) = -1*(cy2(2))*c_h2(i-1)-1*(cy2(3))*c_h2(i-2)+(cx2(1))*c_h2_temp(i)+(cx2(2))*c_h2_temp(i-1)+(cx2(3))*c_h2_temp(i-2);
        end
    end
    for i = 1:length(c)
        c_h2_t(i) = c_h2(i);                        %存回去一開始的暫存器，以便做下一個m的計算
    end
end

figure(5);                                          %Plot figure d in homework(額外畫一張原訊號x[n]圖做效果比對)
subplot(2,1,1);                                                
plot(c);
title('Waveform of signal x[n]');
grid on
xlabel("Sample number (n)");
subplot(2,1,2); 
plot(c_h2);
title('Waveform of signal x[n] after H1 and H2');
grid on
xlabel("Sample number (n)");



%% Plot unit circle start
xp = 1*cos(ang);
yp = 1*sin(ang);
figure(6);              %這張圖是多畫的，畫出unit circle的pole-zero圖形
plot(real(thetas_c),imag(thetas_c),"o",real(thetas_x),imag(thetas_x),"x",xp,yp,'b'); 
axis equal
grid on
title('Pole-zero plot for the filter');
xlabel("Re(z)")
ylabel("Im(z)")

toc;            %Stop timer(計時耗費時間用，Command Window跳出Elapsed time才是程式結束)
                %程式跑的時間約為40秒，根據每台電腦的性能或是angle解析度有所差異(d_ang越小時間越長，反之則時間越短)

%% DTFT function
function [X,x,df]=fftseq(m,ts,df)
% [M,m1,df]=fftseq(m,ts,df)
% [M,m1,df]=fftseq(m,ts) 不輸入df將以最低解析度計算
% m為欲轉換的信號 ts為取樣間隔時間 df為解析度
fs=1/ts;
if nargin == 2
n1=0;
else
n1=fs/df;                                %n1 = 假設頻域上點數
end
n2=length(m);                            %n2 = 實際時域上點數
n=2^(max(nextpow2(n1),nextpow2(n2)));    %取點數最多的那方並找出最靠近2的冪次方的數字
X=fft(m,n);                              %X為FFT轉換後的信號
x=[m,zeros(1,n-n2)];                     %x為將原信號後端補零使取樣個數n為2的冪次方
df=fs/n;                                 %df為此轉換程式真正使用的解析度
end