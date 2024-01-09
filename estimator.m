m1=0; %元件有無noise
m2=0;
m3=0;
m4=0;
m5=0;
Zr2=1000000; %1M
Zr1=1000000; %1M
Zc1=jw*10^-6;%w=?
Zc2=jw*10^-6; 
dT=0.1; %取點間的時間差

X= [N1 N2 N3 N4 V1 Vout]; %上一步的X

A55= -1/((Zr1+m1*N1)*(Zc1+m3*N3))-1/((Zr2+m2*N2)*(Zc1+m3*N3))
A56= 1/((Zr2+m2*N2)*(Zc1+m3*N3));
A65= 1/((Zr2+m2*N2)*(Zc2+m4*N4));
A66= -1/((Zr2+m2*N2)*(Zc2+m4*N4));
A= [1-dT 0 0 0 0 0;
    0 1-dT 0 0 0 0;
    0 0 1-dT 0 0 0;
    0 0 0 1-dT 0 0;
    0 0 0 0 1+dT*A55 A56;
    0 0 0 0 A65 1+dTA66] %𝜌= 1

Z5= dT*u/((Zr1+m1*N1)(Zc1+m3*N3));
Z= transpose([0 0 0 0 Z5 0]);

K56= m5/((Zr1+m1*N1)(Zc1+m3*N3)); %k該是matrix還是一個值？
K= [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 K56;
    0 0 0 0 0 0] %非常不確定
B= transpose([wgn(1, 5, 50) 0]) % power of noise 50dbW, brownian motion process微分是高斯雜訊吧
c=[0 0 0 0 1 0];

Xnew= A*X+ Z+ K56*B;
Y=c*X;