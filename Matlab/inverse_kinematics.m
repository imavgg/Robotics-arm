clear
clc


%================================== test mode =============================
% 
N = [0.0805 0.3571 0.9306];
O = [-0.3572 0.8819 -0.3076];
A = [-0.9306 -0.3076 0.1986];
P = [-7.9883 8.5915 -1.7365];

%================================== cartesian point =======================

nx = N(1);ny = N(2);nz = N(3);                 % define n
ox = O(1);oy = O(2);oz = O(3);                 % define o
ax = A(1);ay = A(2);az = A(3);                 % define a 
px = P(1);py = P(2);pz = P(3);                 % define p
% Cartesian_Point = [nx ox ax px; ny oy ay py; nz oz az pz; 0 0 0 1];

%================================== kinematic table =========================
%four possible answers of [ theta1~theta6]
answer = zeros(4,6); 
d2 = 6.375; 

%================================= solve theta 1 :2sol==========================

theta_1 = zeros(1,4);
FiTheta = zeros(1,4);
FiTheta(1:2) = 180/pi*atan2(d2,((px^2 + py^2 - d2^2)^0.5)); %1
FiTheta(3:4) = 180/pi*atan2(d2,(-(px^2 + py^2 - d2^2)^0.5));%2
for i = 1:1:4
 
    theta_1(i) = 180/pi*atan2(py,px) - FiTheta(i);
end
answer(1:4,1) = theta_1(1:4); 


%================================= solve theta 2  ==========================

theta_2 = zeros(1,4);
for i = 1:1:4 
       m=cos(pi/180*theta_1(i))*px+sin(pi/180*theta_1(i))*py ;
       theta_2(i) = (180/pi)*atan2( m, pz);
end
answer(1:4,2) = theta_2(1:4);

%================================= solve d3==========================
d3=zeros(1,4);
for i = 1:4
d3(i) = cos(pi/180*theta_1(i))*sin(pi/180*theta_2(i))*px+sin(pi/180*theta_1(i))*sin(pi/180*theta_2(i))*py+cos(pi/180*theta_2(i))*pz;
theta3 =0;
end
%================================= solve theta 4 ==========================
theta_4 = zeros(1,4);
for i = 1:2:4 %1,3 input  2 solutions of theta1,theta2
        c1 = cos(pi/180*theta_1(i));
        s1 = sin(pi/180*theta_1(i));
        
        c2 = cos(pi/180*theta_2(i));
        s2 = sin(pi/180*theta_2(i));
                     
  
        m=-s1*ax+c1*ay;
        n=c1*c2*ax+s1*c2*ay-s2*az;

        theta_4(i) = 180/pi*atan2(m,n);
end

for i = 2:2:4 %2,4 input  2 solutions of theta1,theta2
        c1 = cos(pi/180*theta_1(i));
        s1 = sin(pi/180*theta_1(i));
        
        c2 = cos(pi/180*theta_2(i));
        s2 = sin(pi/180*theta_2(i));

        m=-s1*ax+c1*ay;
        n=c1*c2*ax+s1*c2*ay-s2*az;     
        theta_4(i) = 180/pi*atan2(-m,-n);
end
answer(1:4,4) = theta_4(1:4);
%================================= solve theta 6 ==========================

theta_6 = zeros(1,4);
for i = 1:4 %1,3 input  2 solutions of theta1,theta2
        c1 = cos(pi/180*theta_1(i));
        s1 = sin(pi/180*theta_1(i));
        
        c2 = cos(pi/180*theta_2(i));
        s2 = sin(pi/180*theta_2(i));
       
        m = c1*s2*ox+s1*s2*oy+c2*oz;
        n = c1*s2*nx+s1*s2*ny+c2*nz;
        if i == 1 || i==3
            theta_6(i) = 180/pi*atan2(m,-n);

        end

        if i == 4 || i==2
            theta_6(i) = 180/pi*atan2(-m,n);
        end
end

answer(1:4,6) = theta_6(1:4);
%================================= solve theta 5 ==========================

theta_5 = zeros(1,4);
for i = 1:4
        c1 = cos(pi/180*theta_1(i));
        s1 = sin(pi/180*theta_1(i));
        
        c2 = cos(pi/180*theta_2(i));
        s2 = sin(pi/180*theta_2(i));
      
        c4 = cos(pi/180*theta_4(i));
        s4 = sin(pi/180*theta_4(i));
        
        m=ax*(c4*c1*c2-s4*s1)+ay*(c4*c2*s1+c1*s4)-c4*s2*az;
        n=-c1*s2*ax-s1*s2*ay-c2*az;     
            theta_5(i) = 180/pi*atan2(m,-n);

end

answer(1:4,5) = theta_5(1:4);
answer();

%================================= print the answer =======================

 fprintf('  theta1 , theta2 , d3 , theta4  , theta5 , theta6 =',i);


for i = 1:1:4
        fprintf('\n---------------------------------------- \nAnswer(%d)\n',i);   
        fprintf('       %f   %f   %f   %f   %f   %f  , %f  \n\n',theta_1(i), theta_2(i), d3(i), theta_4(i), theta_5(i), theta_6(i));
       
        if (160 <= theta_1(i) ||  theta_1(i) <= -160)
            fprintf(' \n theta 1 is out of range\n'); 
        end

        if (125 <= theta_2(i) ||  theta_2(i) <= -125)
        else
         fprintf('\n theta 2 is out of range\n'); 
        end

        if (30 <= d3(i) ||  d3(i)< -30)
           fprintf('\n d3 is out of range\n'); 
        end

        if (140 <= theta_4(i) ||  theta_4(i) <= -140)
           fprintf('\n theta 4 is out of range\n'); 
        end

        if (100 <= theta_5(i) ||  theta_5(i) <= -100)
           fprintf('\n theta 5 is out of range\n'); 
        end

        if (260 <= theta_6(i) ||  theta_6(i) <= -260)
          fprintf('\n theta 6 is out of range\n'); 
        end   
end
