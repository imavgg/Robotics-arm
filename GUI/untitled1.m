
function varargout = gui_project1(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled1_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

global NOAP;%=zeros(4,4);
global A;%=zeros(1,6);
global B;%=zeros(1,9);

% --- Executes just before untitled1 is made visible.
function untitled1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled1 (see VARARGIN)

% Choose default command line output for untitled1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%  kinematics
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta_1=str2double(get(handles.edit1,'String'));
theta_2=str2double(get(handles.edit3,'String'));
theta_4=str2double(get(handles.edit4,'String'));
theta_5=str2double(get(handles.edit5,'String'));
theta_6=str2double(get(handles.edit6,'String'));
d3=str2double(get(handles.edit2,'String')) ;

theta1=deg2rad(theta_1);
theta2=deg2rad(theta_2);
theta3=0;
theta4=deg2rad(theta_4);
theta5=deg2rad(theta_5);
theta6=deg2rad(theta_6);
%================================= definition ==========================
%definitions of A1~A6 
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

NOAP=T6;
A=[x;y;z;phi(1);theta(1);psi(1)]';
B=[x;y;z;phi(1); theta(1);psi(1);   phi(1)*pi/180 ; theta(1)*pi/180; psi(1)*pi/180 ]';
% 
% guidata(NOAP,T6);
% handles.arrayData=T6;
% guidata(hObject,handles);
set(handles.uitable3,'Data', NOAP );
set(handles.edit7,'String', num2str(B) );

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
% NOAP
% (handles.edit7,'string',mat2str(NOAP));
% set(handles.edit7,'string',mat2str(NOAP));
% --- Executes during object creation, after setting all properties.


function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% set(handles.uitable2,'String',num2str(answer));


% --- Executes when entered data in editable cell(s) in uitable3.
function uitable3_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable3 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


%% inverse kinematics

function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
nx=str2double(get(handles.edit9,'String'));
ny=str2double(get(handles.edit13,'String'));
nz=str2double(get(handles.edit17,'String'));
ox=str2double(get(handles.edit10,'String'));
oy=str2double(get(handles.edit14,'String'));
oz=str2double(get(handles.edit18,'String')) ;
ax=str2double(get(handles.edit11,'String'));
ay=str2double(get(handles.edit15,'String'));
az=str2double(get(handles.edit19,'String'));
px=str2double(get(handles.edit12,'String'));
py=str2double(get(handles.edit16,'String'));
pz=str2double(get(handles.edit20,'String')) ;
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 
% N = [0.0805 0.3571 0.9306];
% O = [-0.3572 0.8819 -0.3076];
% A = [-0.9306 -0.3076 0.1986];
% P = [-7.9883 8.5915 -1.7365];

%================================== cartesian point =======================
% 
% nx = N(1);ny = N(2);nz = N(3);                 % define n
% ox = O(1);oy = O(2);oz = O(3);                 % define o
% ax = A(1);ay = A(2);az = A(3);                 % define a 
% px = P(1);py = P(2);pz = P(3);                 % define p
% % Cartesian_Point = [nx ox ax px; ny oy ay py; nz oz az pz; 0 0 0 1];

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

for i = 1:4
d3 = cos(pi/180*theta_1(i))*sin(pi/180*theta_2(i))*px+sin(pi/180*theta_1(i))*sin(pi/180*theta_2(i))*py+cos(pi/180*theta_2(i))*pz;
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

all2=zeros(4,6);
for i = 1:1:4
       
        if (160 <= theta_1(i) ||  theta_1(i) <= -160)
            fprintf('\n theta 1 is out of range\n');
        end

        if (125 <= theta_2(i) ||  theta_2(i) <= -125)
            fprintf('\n theta 2 is out of range\n');
        end

        if (30 <= d3 ||  d3< -30)
            fprintf('\n d3 is out of range\n');
        end

        if (140 <= theta_4(i) ||  theta_4(i) <= -140)
          fprintf('\n theta 4 is out of range\n'); 

        end

        if (100 <= theta_1(i) ||  theta_1(i) <= -100)
          fprintf('\n theta 5 is out of range\n'); 

        end

        if (260 <= theta_6(i) ||  theta_6(i) <= -260)
          fprintf('\n theta 6 is out of range\n'); 
        end   
        all2(i,:)=[theta_1(i), theta_2(i), d3, theta_4(i), theta_5(i), theta_6(i)];

                
end
set(handles.uitable4,'Data', all2 );




function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable4.
function uitable4_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable4 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
