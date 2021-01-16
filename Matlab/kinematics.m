% Stanford
%  

clear
clc
%================================== test mode =============================
joint_variables=[100,100,100,100,100];  %theta3=0
d3=10;

%================================= input ==========================


theta1=joint_variables(1,1);
theta2=joint_variables(1,2);
theta3 = 0;
theta4=joint_variables(1,3);
theta5=joint_variables(1,4);
theta6=joint_variables(1,5);

theta1=deg2rad(theta1);
theta2=deg2rad(theta2);
theta4=deg2rad(theta4);
theta5=deg2rad(theta5);
theta6=deg2rad(theta6);
%================================= definition ==========================
% definitions of A1~A6 
A1=[cos(theta1) 0 -sin(theta1)      0;
    sin(theta1) 0  cos(theta1)      0;
         0     -1       0           0  ;      
         0      0       0           1         ];
     
A2=[cos(theta2)       0 sin(theta2)   0;
    sin(theta2)       0 -cos(theta2)  0;
         0            1      0        6.375      ;  
         0            0      0        1         ];
     
A3=[cos(theta3) -sin(theta3) 0      0        ;
    sin(theta3)  cos(theta3) 0      0       ; 
         0            0      1      d3       ;
         0            0      0      1         ];
     
     
A4=[cos(theta4) 0 -sin(theta4)      0;
    sin(theta4) 0  cos(theta4)      0;
         0     -1       0           0  ;      
         0      0       0           1         ];
     
A5=[cos(theta5) 0  sin(theta5)      0       ; 
    sin(theta5) 0 -cos(theta5)      0        ;
         0      1       0           0      ;  
         0      0       0           1         ];
     
A6=[cos(theta6) -sin(theta6) 0      0        ;
    sin(theta6)  cos(theta6) 0      0       ; 
         0            0      1      0        ;
         0            0      0      1         ];


T6=A1*A2*A3*A4*A5*A6; %noap
nx=T6(1,1);ny=T6(2,1);nz=T6(3,1);
ox=T6(1,2);oy=T6(2,2);oz=T6(3,2);
ax=T6(1,3);ay=T6(2,3);az=T6(3,3);
px=T6(1,4);py=T6(2,4);pz=T6(3,4);

%================================= solve ==========================


x=px;
y=py;
z=pz;

phi(1)=atan2(ay,ax);
phi(2)=atan2(-ay,-ax);


theta(1)=atan2( (cos(phi(1))*ax+sin(phi(1))*ay) ,az);
theta(2)=atan2( (cos(phi(2))*ax+sin(phi(2))*ay)  ,az);


psi(1)=atan2( (-sin(phi(1))*nx+cos(phi(1))*ny) , (-sin(phi(1))*ox+cos(phi(1))*oy) );
psi(2)=atan2( (-sin(phi(2))*nx+cos(phi(2))*ny) , (-sin(phi(2))*ox+cos(phi(2))*oy) );


phi=[phi(1) phi(2)]*180/pi;
theta=[theta(1) theta(2)]*180/pi;
psi=[psi(1) psi(2)]*180/pi;
%================================= print the answer =======================

disp("(n,o,a,p)=");
disp(T6);
disp("-------------------------------------------");
disp("(x,y,z,phi,theta,psi)=");
disp([x;y;z;phi(1);theta(1);psi(1)]');

disp("-------------------------------------------");
disp("cartesian point=");
disp([x;y;z;phi(1); theta(1);psi(1);   phi(1)*pi/180 ; theta(1)*pi/180; psi(1)*pi/180 ]');

