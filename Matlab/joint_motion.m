function output = joint_motion(A,B,C)
%these can be changed
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
joint_a=inverse_kinematics(A_n,A_o,A_a,A_p);
 % use Project1 : Inverse kinematics :get the same solution in report
joint_a=joint_a(2,:);
B_n=[B(1:3,1)];
B_o=[B(1:3,2)];
B_a=[B(1:3,3)];
B_p=[B(1:3,4)];
joint_b=inverse_kinematics(B_n,B_o,B_a,B_p);
 % use Project1 : Inverse kinematics :get the same solution in report
joint_b=joint_b(2,:);

C_n=[C(1:3,1)];  
C_o=[C(1:3,2)];
C_a=[C(1:3,3)];
C_p=[C(1:3,4)];
joint_c=inverse_kinematics(C_n,C_o,C_a,C_p);
 % use Project1 : Inverse kinematics :get the same solution in report
joint_c=joint_c(4,:);

theta_A_a=joint_a';                                                            
theta_B_a=joint_b';                              
theta_C_a=joint_c';                       


% AA'
thetA_a2=theta_A_a+(theta_B_a-theta_A_a)/0.5*(0.5-0.2);              
d_B=theta_A_a-theta_B_a;
d_C=theta_C_a-theta_B_a;
count=1;
step=0.002;
time_start=-0.5;
time_end = -0.2;
for t=time_start:step:time_end 
    dtheta_A_a(:,count)=theta_B_a-d_B/0.5*t;   %角度     
    domegA_a(:,count)=-d_B/0.5;                %角速度
    dalphA_a(:,count)=[0;0;0;0;0;0];         %線性區 角加速度為0
    count=count+1;
end

count=1;
step=0.002;
time_start=-0.198;
time_end =0.198;
for t=time_start:step:time_end
    h=(t+0.2)/0.4;
    dtheta_B_a(:,count)=(0.2/0.5)*((d_C+d_B)*(2-h)*(h.^2)-2*d_B)*h+d_B*(0.2/0.5)+theta_B_a;      %角度     
    omega_B(:,count)=((d_C+d_B)*(1.5-h)*2*(h.^2)-d_B)/0.5;                                    %角速度
    alphaB_a(:,count)=(d_C+d_B)*(1-h)*3*h/(0.2*0.5);                                      %非線性區 角加速度不為0                   
    count=count+1;
end


% C'C
count=1;
step=0.002;
time_start=0.2;
time_end = 0.5;
for t=time_start:step:time_end
    h=t/0.5;
    dtheta_C_a(:,count)=d_C*h+theta_B_a;                %角度                          
    domegC_a(:,count)=d_C/0.5;                           %角速度
    dalphC_a(:,count)=[0;0;0;0;0;0];                    %線性區 角加速度為0
    count=count+1;
end

% AC joint
figure(1)
step=0.002;
time_start=-0.5;
time_end = 0.5;
t=time_start:step:time_end;
theta1=[dtheta_A_a(1,:) dtheta_B_a(1,:) dtheta_C_a(1,:)];                 
subplot(3,2,1);
plot(t,theta1);
title('Joint value (degree or cm)');
ylabel('joint1');
theta2=[dtheta_A_a(2,:) dtheta_B_a(2,:) dtheta_C_a(2,:)];                    
subplot(3,2,2);
plot(t,theta2);
ylabel('joint2');
theta3=[dtheta_A_a(3,:) dtheta_B_a(3,:) dtheta_C_a(3,:)];                    
subplot(3,2,3);
plot(t,theta3);
ylabel('joint3');
theta4=[dtheta_A_a(4,:) dtheta_B_a(4,:) dtheta_C_a(4,:)];                   
subplot(3,2,4);
plot(t,theta4);
ylabel('joint4');
theta5=[dtheta_A_a(5,:) dtheta_B_a(5,:) dtheta_C_a(5,:)];                   
subplot(3,2,5);
plot(t,theta5);
ylabel('joint5');
theta6=[dtheta_A_a(6,:) dtheta_B_a(6,:) dtheta_C_a(6,:)];                     
subplot(3,2,6);
plot(t,theta6);
ylabel('joint6');
% 
% % AC velocity
% 
figure(2)
subplot(3,2,1);
plot(t,[domegA_a(1,:) omega_B(1,:) domegC_a(1,:)]);                 
title('Velocity (degree/s or cm/s)');
ylabel('joint1');

subplot(3,2,2);
plot(t,[domegA_a(2,:) omega_B(2,:) domegC_a(2,:)]);    
ylabel('joint2');

subplot(3,2,3);
plot(t,[domegA_a(3,:) omega_B(3,:) domegC_a(3,:)]);     
ylabel('joint3');

subplot(3,2,4);
plot(t,[domegA_a(4,:) omega_B(4,:) domegC_a(4,:)]);      
ylabel('joint4');

subplot(3,2,5);
plot(t,[domegA_a(5,:) omega_B(5,:) domegC_a(5,:)]);       
ylabel('joint5');

subplot(3,2,6);
plot(t,[domegA_a(6,:) omega_B(6,:) domegC_a(6,:)]);          
ylabel('joint6');
% 
% 
% 
% 
% % AC acceleration
figure(3)
subplot(3,2,1);
ylabel('joint1');
plot(t,[dalphA_a(1,:) alphaB_a(1,:) dalphC_a(1,:)]);
title('Acceleration (degree/s^2 or cm/s^2)'); 

subplot(3,2,2);
plot(t,[dalphA_a(2,:) alphaB_a(2,:) dalphC_a(2,:)]);
ylabel('joint2');


subplot(3,2,3);
plot(t,[dalphA_a(3,:) alphaB_a(3,:) dalphC_a(3,:)]);
ylabel('joint3');

subplot(3,2,4);
plot(t,[dalphA_a(4,:) alphaB_a(4,:) dalphC_a(4,:)]);
ylabel('joint4');

subplot(3,2,5);
plot(t,[dalphA_a(5,:) alphaB_a(5,:) dalphC_a(5,:)]);
ylabel('joint5');

subplot(3,2,6);
plot(t,[dalphA_a(6,:) alphaB_a(6,:) dalphC_a(6,:)]);
ylabel('joint6');

output= zeros(6,1); %initialize

% AA'
step=0.002;
time_start=-0.5;
time_end = -0.2;
count=1;
for t=time_start:step:time_end
    % get [ px py pz ]from project1: forward_kinematics  
     AA = forward_kinematics(dtheta_A_a(:,count));     
     x1(count)=AA(1);
     y1(count)=AA(2);
     z1(count)=AA(3);
     output = dtheta_A_a(:,count);
    count=count+1;
end

% A'C'
step=0.002;
time_start=-0.198;
time_end = 0.198;
count=1;
for t=time_start:step:time_end
    % get [ px py pz ]from project1: forward_kinematics  
     AC = forward_kinematics(dtheta_B_a(:,count));
     x2(count)=AC(1);
     y2(count)=AC(2);
     z2(count)=AC(3);
    output = dtheta_B_a(:,count);
    count=count+1;
end

% CC'
step=0.002;
time_start=0.2;
time_end = 0.5;
count=1;
for t=time_start:step:time_end
% get [ px py pz ]from project1: forward_kinematics  
   CC =forward_kinematics(dtheta_C_a(:,count));
    x3(count)=CC(1);
    y3(count)=CC(2);
    z3(count)=CC(3);
    output =dtheta_C_a(:,count);
    count=count+1;
end


% plot
figure(4)
plot3(x1,y1,z1,x2,y2,z2,x3,y3,z3);
xlabel('x(cm)');
ylabel('y(cm)');
zlabel('z(cm)');
text(20,10,-10,'A(20,10,-10)');
text(20,-5,10,'B(20,-5,10)');
text(-10,15,25,'C(-10,15,25)');
grid
title('3D path of Joint Motion')
%PLOT AXIS
hold on;
% A
plot3([A_p(1),A_p(1)],[A_p(2),A_p(2)],[A_p(3),A_p(3)-A_n(3)],'r');
hold on;
plot3([A_p(1),A_p(1)+A_o(1)*5],[A_p(2),A_p(2)+A_o(2)*5],[A_p(3),A_p(3)+A_o(3)*5],'g');
hold on;
plot3([A_p(1),A_p(1)+A_a(1)*5],[A_p(2),A_p(2)+A_a(2)*5],[A_p(3),A_p(3)+A_a(3)*5],'b');
hold on;

% B
plot3([B_p(1),B_p(1)+B_n(1)*5],[B_p(2),B_p(2)+B_n(2)*5],[B_p(3),B_p(3)+B_n(3)*5],'r');
hold on;
plot3([B_p(1),B_p(1)+B_o(1)*5],[B_p(2),B_p(2)+B_o(2)*5],[B_p(3),B_p(3)+B_o(3)*5],'g');
hold on;
plot3([B_p(1),B_p(1)+B_a(1)*5],[B_p(2),B_p(2)+B_a(2)*5],[B_p(3),B_p(3)+B_a(3)*5],'b');
hold on;

% C
plot3([C_p(1),C_p(1)+C_n(1)*5],[C_p(2),C_p(2)+C_n(2)*5],[C_p(3),C_p(3)+C_n(3)*5],'r');
hold on;
plot3([C_p(1),C_p(1)+C_o(1)*5],[C_p(2),C_p(2)+C_o(2)*5],[C_p(3),C_p(3)+C_o(3)*5],'g');
hold on;
plot3([C_p(1),C_p(1)+C_a(1)*5],[C_p(2),C_p(2)+C_a(2)*5],[C_p(3),C_p(3)+C_a(3)*5],'b');
% hold on;