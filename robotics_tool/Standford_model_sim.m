
jv= input('joint varialble =\n');
theta1=jv(1);theta2=jv(2);theta4=jv(4);theta5=jv(5);theta6=jv(6);d3=jv(3);
%=================== test mode=====================================
%input kinematics

% theta1=100;theta2=100;theta4=100;theta5=100;theta6=100;d3=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=Link(zeros(6,5));
% bot=SerialLink([L1 L2 L3 L4 L5 L6),'name','myrobot')
%             th    d       a    alpha
L(1) = Link([ 0    0       0   -pi/2 ]);
L(2) = Link([ 0    6.375 0    pi/2   ]);
L(3) = Link([ 0 0  10    0     0     ]);  % PRISMATIC link
L(4) = Link([ 0     0    0   -pi/2   ]);
L(5) = Link([ 0     0    0    pi/2   ]);
L(6) = Link([ 0     0    0    0      ]);

%These are dynamcis of model
L(3).qlim = [-30 30] * 0.0254;%inch to meter

L(1).qlim = [-160 160]*pi/180;
L(2).qlim = [-125 125]*pi/180;
L(4).qlim = [-140 140]*pi/180;
L(5).qlim = [-100 100]*pi/180;
L(6).qlim = [-260 260]*pi/180;


L(1).m = 9.29;
L(2).m = 5.01;
L(3).m = 4.25;
L(4).m = 1.08;
L(5).m = 0.630;
L(6).m = 0.51;

L(1).Jm = 0.953;
L(2).Jm = 2.193;
L(3).Jm = 0.782;
L(4).Jm = 0.106;
L(5).Jm = 0.097;
L(6).Jm = 0.020;

L(1).G = 1;
L(2).G = 1;
L(3).G = 1;
L(4).G = 1;
L(5).G = 1;
L(6).G = 1;

L(1).I = [0.276   0.255   0.071   0   0   0];
L(2).I = [0.108   0.018   0.100   0   0   0];
L(3).I = [2.51    2.51    0.006   0   0   0 ];
L(4).I = [0.002   0.001   0.001   0   0   0 ];
L(5).I = [0.003   0.0004  0       0   0   0];
L(6).I = [0.013   0.013   0.0003  0   0   0];

L(1).r = [0    0.0175 -0.1105];
L(2).r = [0   -1.054  0];
L(3).r = [0    0      -6.447];
L(4).r = [0    0.092  -0.054];
L(5).r = [0    0.566   0.003];
L(6).r = [0    0       1.554];

% output 
stanf = SerialLink(L, 'name', 'Stanford arm');
stanf.plot([theta_1, theta_2, d3, theta_4, theta_5, theta_6],'floorlevel',20);


