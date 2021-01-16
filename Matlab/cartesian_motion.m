
clc;
clear;
A=[ 0 1 0 20;
        0 0 -1 10;
       -1 0 0 -10;
        0 0 0 1];
B=[-1 0 0 20;
        0 -1 0 -5;
        0 0  1 10;
        0 0 0 1];
C=[ 1 0 0 -10;
        0 0 -1 15;
        0 1 0 25;                                                                                            
        0 0 0 1];
A_n=[A(1:3,1)];
A_o=[A(1:3,2)];
A_a=[A(1:3,3)];
A_p=[A(1:3,4)];

B_n=[B(1:3,1)];
B_o=[B(1:3,2)];
B_a=[B(1:3,3)];
B_p=[B(1:3,4)];

C_n=[C(1:3,1)];
C_o=[C(1:3,2)];
C_a=[C(1:3,3)];
C_p=[C(1:3,4)];
%
% 計算線性區域 A到A'點 的 x y z psi theta phi

x=A_n'*(B_p-A_p);
y=A_o'*(B_p-A_p);
z=A_a'*(B_p-A_p);
psi=atan2(A_o'*B_a,A_n'*B_a);
theta=atan2(sqrt((A_n'*B_a)^2+(A_o'*B_a)^2),A_a'*B_a);
V = 1-cos(theta);
Sin_phi=-sin(psi)*cos(psi)*(V)*(A_n'*B_n)+((cos(psi))^2*(V)+cos(theta))*(A_o'*B_n)-sin(psi)*sin(theta)*(A_a'*B_n);
Cos_phi=-sin(psi)*cos(psi)*(V)*(A_n'*B_o)+((cos(psi))^2*(V)+cos(theta))*(A_o'*B_o)-sin(psi)*sin(theta)*(A_a'*B_o);
phi=atan2(Sin_phi,Cos_phi);
% 計算出從A到A'中的每個點
count=1;
step=0.002;
time_start=0.0;
time_end = 0.3;

for t=time_start:step:time_end
%     h=(t+0.5)/0.5;  %%%%%%%
     h=t/0.5;
    d_x=x*h;
    d_y=y*h;
    d_z=z*h;
    d_psi=psi;
    d_theta=theta*h;
    d_phi=phi*h;   
    
    Sin_psi=sin(psi);
    Cos_psi=cos(psi);
    Sin_theta=sin(d_theta);
    Cos_theta=cos(d_theta);
    V_theta=1-Cos_theta;
    Sin_phi=sin(d_phi);
    Cos_phi=cos(d_phi);
    
    T_r=[   1 0 0 d_x;
                0 1 0 d_y;
                0 0 1 d_z;
                0 0 0 1];

    R_ar=[Sin_psi^2*V_theta+Cos_theta  -Sin_psi*Cos_psi*V_theta    Cos_psi*Sin_theta  0;
               -Sin_psi*Cos_psi*V_theta    Cos_psi^2*V_theta+Cos_theta  Sin_psi*Sin_theta   0;
                 -Cos_psi*Sin_theta            -Sin_psi*Sin_theta           Cos_theta      0;
                      0                         0                0               1];
    R_or=[Cos_phi  -Sin_phi      0  0;
                 Sin_phi   Cos_phi      0  0;
                    0       0       1  0;
                    0       0       0  1];
        
    Dr=T_r*R_ar*R_or; 
 
    A_p_B(:,:,count)=A*Dr;
    xA_B(:,count)=A_p_B(1,4,count);
    yA_B(:,count)=A_p_B(2,4,count);
    zA_B(:,count)=A_p_B(3,4,count);  
    count=count+1;
end

% 計算非線性區域 A'到C'點 的 x y z psi theta phi
A2=A_p_B(:,:,count-1);        % 線性區端  (A-B) 最後一點為非線性區段(A'-C')的起始點
A_n2=[A2(1,1);A2(2,1);A2(3,1)];
A_o2=[A2(1,2);A2(2,2);A2(3,2)];
A_a2=[A2(1,3);A2(2,3);A2(3,3)];
A_p2=[A2(1,4);A2(2,4);A2(3,4)];

xA=B_n'*(A_p2-B_p);
yA=B_o'*(A_p2-B_p);
zA=B_a'*(A_p2-B_p);
psiA=atan2(B_o'*A_a2,B_n'*A_a2);
thetA_a=atan2(sqrt((B_n'*A_a2)^2+(B_o'*A_a2)^2),B_a'*A_a2);
V=1-cos(thetA_a);
SinphiA=-sin(psiA)*cos(psiA)*(V)*(B_n'*A_n2)+((cos(psiA))^2*(V)+cos(thetA_a))*(B_o'*A_n2)-sin(psiA)*sin(thetA_a)*(B_a'*A_n2);
CosphiA=-sin(psiA)*cos(psiA)*(V)*(B_n'*A_o2)+((cos(psiA))^2*(V)+cos(thetA_a))*(B_o'*A_o2)-sin(psiA)*sin(thetA_a)*(B_a'*A_o2);
phiA=atan2(SinphiA,CosphiA);

% 計算非線性區域 B到C點 的 x y z psi theta phi
xC=B_n'*(C_p-B_p);
yC=B_o'*(C_p-B_p);
zC=B_a'*(C_p-B_p);
psiC=atan2(B_o'*C_a,B_n'*C_a);
thetC_a=atan2(sqrt((B_n'*C_a)^2+(B_o'*C_a)^2),B_a'*C_a);
V= 1-cos(thetC_a);
SinphiC=-sin(psiC)*cos(psiC)*(V)*(B_n'*C_n)+((cos(psiC))^2*(V)+cos(thetC_a))*(B_o'*C_n)-sin(psiC)*sin(thetC_a)*(B_a'*C_n);
CosphiC=-sin(psiC)*cos(psiC)*(V)*(B_n'*C_o)+((cos(psiC))^2*(V)+cos(thetC_a))*(B_o'*C_o)-sin(psiC)*sin(thetC_a)*(B_a'*C_o);
phiC=atan2(SinphiC,CosphiC);

if abs(psiC-psiA)>pi/2
    psiA=psiA+pi;
    thetA_a=-thetA_a;
end
% 計算出從A'到C'中的每個點
count=1;
step = 0.002;
time_start=0.302;
time_end=0.698;
for t=time_start:0.002:time_end
    h=(t-0.30)/0.4;  %%%%%%%
    dx_B=((xC*0.2/0.5+xA)*(2-h)*h^2-2*xA)*h+xA;
    dy_B=((yC*0.2/0.5+yA)*(2-h)*h^2-2*yA)*h+yA;
    dz_B=((zC*0.2/0.5+zA)*(2-h)*h^2-2*zA)*h+zA;
    dpsi_B=(psiC-psiA)*h+psiA;
    dtheta_B=((thetC_a*0.2/0.5+thetA_a)*(2-h)*h^2-2*thetA_a)*h+thetA_a;
    dphi_B=((phiC*0.2/0.5+phiA)*(2-h)*h^2-2*phiA)*h+phiA;
  
    Sin_psi_B=sin(dpsi_B);
    Cos_psi_B=cos(dpsi_B);
    Sin_theta_B=sin(dtheta_B);
    Cos_theta_B=cos(dtheta_B);
    V_theta_B=1-Cos_theta_B;
    Sin_phi_B=sin(dphi_B);
    Cos_phi_B=cos(dphi_B);
    
  
    T_r_B=[1 0 0 dx_B;
          0 1 0 dy_B;
          0 0 1 dz_B;
          0 0 0   1 ];
 
    R_ar_B=[Sin_psi_B^2*V_theta_B+Cos_theta_B  -Sin_psi_B*Cos_psi_B*V_theta_B         Cos_psi_B*Sin_theta_B   0;
         -Sin_psi_B*Cos_psi_B*V_theta_B        Cos_psi_B^2*V_theta_B+Cos_theta_B     Sin_psi_B*Sin_theta_B   0;
         -Cos_psi_B*Sin_theta_B               -Sin_psi_B*Sin_theta_B                     Cos_theta_B         0;
              0                                              0                                0              1];
          
    R_or_B=[Cos_phi_B  -Sin_phi_B      0  0;
         Sin_phi_B   Cos_phi_B      0  0;
            0       0           1  0;
            0       0           0  1];
        
    Dr_B=T_r_B*R_ar_B*R_or_B;   
    
    p_B(:,:,count)=B*Dr_B;
    x_B(:,count)=p_B(1,4,count);
    y_B(:,count)=p_B(2,4,count);
    z_B(:,count)=p_B(3,4,count);  
    count=count+1;
end

% 計算出從C'到C中的每個點
count=1;
for t=0.7:0.002:1.0
    h=(t-0.7)/0.3;   %%%%%%%
    dx_C=xC*h;
    dy_C=yC*h;
    dz_C=zC*h;
    dpsi_C=psiC;
    dtheta_C=thetC_a*h;
    dphi_C=phiC*h;
    %=====================================
    Sin_psi_C=sin(dpsi_C);
    Cos_psi_C=cos(dpsi_C);
    Sin_theta_C=sin(dtheta_C);
    Cos_theta_C=cos(dtheta_C);
    V_theta_C=1-Cos_theta_C;
    Sin_phi_C=sin(dphi_C);
    Cos_phi_C=cos(dphi_C);
    
    T_r_C=[1 0 0 dx_C;
                  0 1 0 dy_C;
                  0 0 1 dz_C;
                  0 0 0   1 ];
    
    R_ar_C=[Sin_psi_C^2*V_theta_C+Cos_theta_C  -Sin_psi_C*Cos_psi_C*V_theta_C        Cos_psi_C*Sin_theta_C   0;
                    -Sin_psi_C*Cos_psi_C*V_theta_C        Cos_psi_C^2*V_theta_C+Cos_theta_C     Sin_psi_C*Sin_theta_C   0;
                     -Cos_psi_C*Sin_theta_C               -Sin_psi_C*Sin_theta_C                     Cos_theta_C      0;
                          0                                          0                      0               1];
    R_or_C=[Cos_phi_C  -Sin_phi_C      0  0;
                     Sin_phi_C   Cos_phi_C      0  0;
                        0       0           1  0;
                        0       0           0  1];
    Dr_C=T_r_C*R_ar_C*R_or_C;   
 
    p_C(:,:,count)=B*Dr_C;
    x_C(:,count)=p_C(1,4,count);
    y_C(:,count)=p_C(2,4,count);
    z_C(:,count)=p_C(3,4,count);  
    count=count+1;
end

% PLOT
% POS

X=[xA_B x_B x_C];
Y=[yA_B y_B y_C];
Z=[zA_B z_B z_C];
time_start = 0.0;
time_end   = 1.0;
step=0.002;
t=time_start:step:time_end;
figure(1)
subplot(3,1,1);
plot(t,X);
title('position of x');
grid
subplot(3,1,2);
plot(t,Y);
title('position of y');
grid
subplot(3,1,3);
plot(t,Z);
title('position of z');
grid

time_step=0.002;
% VEL
d_t=t(2:501);
d_X=diff(X)/time_step;
d_Y=diff(Y)/time_step;
d_Z=diff(Z)/time_step;
figure(2)
subplot(3,1,1);
plot(d_t,d_X);
title('velocity of x');
grid
subplot(3,1,2);
plot(d_t,d_Y);
title('velocity of y');
grid
subplot(3,1,3);
plot(d_t,d_Z);
title('velocity of z');
grid


% ACC

dd_t=t(3:501);
dd_X=diff(d_X)/time_step;   
dd_Y=diff(d_Y)/time_step;
dd_Z=diff(d_Z)/time_step;
figure(3)
subplot(3,1,1);
plot(dd_t,dd_X);
title('acceleration of x');
grid
subplot(3,1,2);
plot(dd_t,dd_Y);
title('acceleration of y');
grid
subplot(3,1,3);
plot(dd_t,dd_Z);
title('acceleration of z');
grid

% PLOT
figure(4)
plot3(xA_B,yA_B,zA_B,x_B,y_B,z_B,x_C,y_C,z_C);
xlabel('x(cm)');
ylabel('y(cm)');
zlabel('z(cm)');
text(20,10,-10,'A(20,10,-10)');
text(5,-15,10,'B(5,-15,10)');
text(-10,15,25,'C(-10,15,25)');

grid;
title('3D path of Cartesion Motion');


%plot AXIS
hold on;
plot3([A_p(1),A_p(1)+A_n(1)*5],[A_p(2),A_p(2)+A_n(2)*5],[A_p(3),A_p(3)+A_n(3)*5],'r');
hold on;
plot3([A_p(1),A_p(1)+A_o(1)*5],[A_p(2),A_p(2)+A_o(2)*5],[A_p(3),A_p(3)+A_o(3)*5],'g');
hold on;
plot3([A_p(1),A_p(1)+A_a(1)*5],[A_p(2),A_p(2)+A_a(2)*5],[A_p(3),A_p(3)+A_a(3)*5],'b');
hold on;
plot3([B_p(1),B_p(1)+B_n(1)*5],[B_p(2),B_p(2)+B_n(2)*5],[B_p(3),B_p(3)+B_n(3)*5],'r');
hold on;
plot3([B_p(1),B_p(1)+B_o(1)*5],[B_p(2),B_p(2)+B_o(2)*5],[B_p(3),B_p(3)+B_o(3)*5],'g');
hold on;
plot3([B_p(1),B_p(1)+B_a(1)*5],[B_p(2),B_p(2)+B_a(2)*5],[B_p(3),B_p(3)+B_a(3)*5],'b');
hold on;
plot3([C_p(1),C_p(1)+C_n(1)*5],[C_p(2),C_p(2)+C_n(2)*5],[C_p(3),C_p(3)+C_n(3)*5],'r');
hold on;
plot3([C_p(1),C_p(1)+C_o(1)*5],[C_p(2),C_p(2)+C_o(2)*5],[C_p(3),C_p(3)+C_o(3)*5],'g');
hold on;
plot3([C_p(1),C_p(1)+C_a(1)*5],[C_p(2),C_p(2)+C_a(2)*5],[C_p(3),C_p(3)+C_a(3)*5],'b');
hold on;

%AA'
s=1;
for j1=0:0.002:0.3
    hold on;
    plot3([xA_B(s),xA_B(s)+A_p_B(1,3,s)/20],[yA_B(s),yA_B(s)+A_p_B(2,3,s)/20],[zA_B(s),zA_B(s)+A_p_B(3,3,s)/20]);
    s=s+1;
end

%A'B

s=1;                                          
for j2=0.302:0.002:0.698
    hold on;
    plot3([x_B(s),x_B(s)+p_B(1,3,s)/20],[y_B(s),y_B(s)+p_B(2,3,s)/20],[z_B(s),z_B(s)+p_B(3,3,s)/20]);
    s=s+1;
end

%BC
s=1;      
for j3=0.70:0.002:1.00
    hold on;
    plot3([x_C(s),x_C(s)+p_C(1,3,s)/20],[y_C(s),y_C(s)+p_C(2,3,s)/20],[z_C(s),z_C(s)+p_C(3,3,s)/20]);
    s=s+1;
end
% 在路徑圖中標上方向以及ABC三點的 座標軸
