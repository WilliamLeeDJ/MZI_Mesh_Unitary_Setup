%Example
%Matrix_Data_Refer_From:
%        F. Shokraneh, M. S. Nezami and O. Liboiron-Ladouceur, 
%        "Theoretical and Experimental Analysis of a 4 × 4 Reconfigurable MZI-Based Linear Optical Processor," 
%        in Journal of Lightwave Technology, vol. 38, no. 6, pp. 1258-1267, 15 March15, 2020, doi: 10.1109/JLT.2020.2966949.
%clear all;    
% APP_TSU4 =...
%         [-0.2341+0.0030j -0.1011+0.1765j 0.5216+0.4673j 0.1664+0.6210j;...
%          0.0953+0.2949j 0.7782+0.2674j 0.1030+0.3555j -0.2120-0.2120j;...
%          0.7987+0.2694j -0.1064+0.2709j -0.1671-0.0954j 0.0357+0.4080j;...
%          0.2852+0.2393j -0.3399-0.2852j 0.1006+0.5704j 0.3290-0.4698j];   
%  [Phi,Theta] = APP_TSU4_To_Angle_Calculate(APP_TSU4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi,Theta] = APP_TSU4_To_Angle_Calculate(APP_TSU4)
%%%-----------Adjustment-------------
    Building_Parameter_Use_Num = 2;
    Building_Parameter_Use_Num_2 = 1;
%%%----------------------------------   
    Venify_APP_TSU4 = APP_TSU4; 
    APP_TSU4_Dagger = transpose(conj(APP_TSU4));
    APP_SU4_Unitary = APP_TSU4*APP_TSU4_Dagger; 
    Check_Special = round(det(APP_SU4_Unitary));
    Theta_Sample = 0:0.01:2*pi;
    N = length(Theta_Sample);
    Phi = zeros(1,6);
    Theta = zeros(1,6);
if(Check_Special ~= 1)
    disp('    -------ATTENTION-------');
    disp('Matrix is not a Special Unitary Matrix');
    disp('-------Calculation Abort-------')
else
    disp('    -------ATTENTION-------');
    disp('Matrix is Special Unitary Matrix');
    disp('-------Calculation Continue-------');
    Equation_Setting_Matrix = APP_TSU4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  THETA1      
    MZI_Theta_1 = abs(2*atan((-1*APP_TSU4(4,2))/(APP_TSU4(4,1))));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [DMZI_R] = DMZI_Setup(MZI_Theta_1,pi/3);
    DMZI_Replace = transpose(conj(DMZI_R));
    DMZI_H44 = [DMZI_Replace(1,1) DMZI_Replace(1,2) 0 0;DMZI_Replace(2,1) DMZI_Replace(2,2) 0 0;0 0 1 0;0 0 0 1];
    Middle_APP = APP_TSU4*DMZI_H44;
    APP_TSU4(4,1) = Middle_APP(4,1);
    APP_TSU4(1,2) = Middle_APP(1,2);
    APP_TSU4(2,2) = Middle_APP(2,2);
    APP_TSU4(3,2) = Middle_APP(3,2);
    APP_TSU4(4,2) = Middle_APP(4,2);
    Equation_Setting_Matrix_MZI1_After = APP_TSU4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THETA2
    Theta_Judge2 = zeros(1,N);
    for i = 1:N
     Theta_Judge2(i) = abs(APP_TSU4(4,2)*sin(Theta_Sample(i)/2)+APP_TSU4(4,3)*cos(Theta_Sample(i)/2));
    end
    Evluate2 = min(Theta_Judge2);
    Position = abs(Theta_Judge2)==Evluate2;
    MZI_Theta_2 = Theta_Sample(Position);
    [DMZI_R] = DMZI_Setup(MZI_Theta_2,pi/3);
    DMZI_Replace = transpose(conj(DMZI_R));
    DMZI_H44 = [1 0 0 0;0 DMZI_Replace(1,1) DMZI_Replace(1,2) 0;0 DMZI_Replace(2,1) DMZI_Replace(2,2) 0;0 0 0 1];
    Middle_APP2 = APP_TSU4*DMZI_H44;
    APP_TSU4(4,2) = Middle_APP2(4,2);
    APP_TSU4(1,3) = Middle_APP2(1,3);
    APP_TSU4(2,3) = Middle_APP2(2,3);
    APP_TSU4(3,3) = Middle_APP2(3,3);
    APP_TSU4(4,3) = Middle_APP2(4,3);     
    Equation_Setting_Matrix_MZI2_After = APP_TSU4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  THETA4  
    Theta_Judge4 = zeros(1,N);
    for i = 1:N
     Theta_Judge4(i) = abs(APP_TSU4(4,3)*sin(Theta_Sample(i)/2)+APP_TSU4(4,4)*cos(Theta_Sample(i)/2));
    end
    Evluate4 = min(Theta_Judge4);
    Position = find(abs(Theta_Judge4)==Evluate4);
    MZI_Theta_4 = Theta_Sample(Position);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Equations
    %%%%%%%%%%%%%%%%%%%%%%%%%% MZI 1 %%%%%%%%%%%%%%%%%%%%%%%%%%
    MZI1_H_Ahead = -1j*exp(-1j*(MZI_Theta_1/2));
    MZI1_H11_sine = MZI1_H_Ahead*sin(MZI_Theta_1/2);
    MZI1_H21_cosine = MZI1_H_Ahead*cos(MZI_Theta_1/2);
    MZI1_H12 = MZI1_H_Ahead*cos(MZI_Theta_1/2);
    MZI1_H22 = MZI1_H_Ahead*(-sin(MZI_Theta_1/2));
    %%%%%%%%%%%%%%%%%%%%%%%%%% MZI 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
    MZI2_H_Ahead = -1j*exp(-1j*(MZI_Theta_2/2));
    MZI2_H11_sine = MZI2_H_Ahead*sin(MZI_Theta_2/2);
    MZI2_H21_cosine = MZI2_H_Ahead*cos(MZI_Theta_2/2);
    MZI2_H12 = MZI2_H_Ahead*cos(MZI_Theta_2/2);
    MZI2_H22 = MZI2_H_Ahead*(-sin(MZI_Theta_2/2));    
    %%%%%%%%%%%%%%%%%%%%%%%%%% MZI 4 %%%%%%%%%%%%%%%%%%%%%%%%%%
    MZI4_H_Ahead = -1j*exp(-1j*(MZI_Theta_4/2));
    MZI4_H11_sine = MZI4_H_Ahead*sin(MZI_Theta_4/2);
    MZI4_H21_cosine = MZI4_H_Ahead*cos(MZI_Theta_4/2);
    MZI4_H12 = MZI4_H_Ahead*cos(MZI_Theta_4/2);
    MZI4_H22 = MZI4_H_Ahead*(-sin(MZI_Theta_4/2));      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %exp(-j phi1)   X1
    %exp(-j phi2)   X2
    %exp(-j phi4)   X4
    %exp(-j Theta3) X3
    %exp(-j Theta5) X5 
%     syms X1 X2 X4 X3 X5
    X1 = sym('X1');
    X2 = sym('X2');
    X3 = sym('X3');
    X4 = sym('X4');
    X5 = sym('X5');
    
    EQ_DMZI1 = [X1*MZI1_H11_sine MZI1_H12;...
             X1*MZI1_H21_cosine MZI1_H22];
    EQ_DMZI1_H44 = [EQ_DMZI1(1,1) EQ_DMZI1(1,2) 0 0;EQ_DMZI1(2,1) EQ_DMZI1(2,2) 0 0;0 0 1 0;0 0 0 1];    
    %-------------GET Cofficienet--------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EQ_DMZI2 = [X2*MZI2_H11_sine MZI2_H12;...
             X2*MZI2_H21_cosine MZI2_H22];
    EQ_DMZI2_H44 = [1 0 0 0;0 EQ_DMZI2(1,1) EQ_DMZI2(1,2) 0;0 EQ_DMZI2(2,1) EQ_DMZI2(2,2) 0;0 0 0 1];
    %-------------GET Cofficienet--------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    EQ_DMZI4 = [X4*MZI4_H11_sine MZI4_H12;...
             X4*MZI4_H21_cosine MZI4_H22];    
    EQ_DMZI4_H44 = [1 0 0 0;0 1 0 0;0 0 EQ_DMZI4(1,1) EQ_DMZI4(1,2);0 0 EQ_DMZI4(2,1) EQ_DMZI4(2,2)];       

    %-------------GET Cofficienet--------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %first_Step MZI equation settings
    EQ_MZI1_MIDDLE = Equation_Setting_Matrix*EQ_DMZI1_H44;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MZI_Equations(1,1) = EQ_MZI1_MIDDLE(1,1);
    MZI_Equations(2,1) = EQ_MZI1_MIDDLE(2,1);
    MZI_Equations(3,1) = EQ_MZI1_MIDDLE(3,1);
    EQ_MZI2_MIDDLE = Equation_Setting_Matrix_MZI1_After*EQ_DMZI2_H44;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    MZI_Equations(1,2) = EQ_MZI2_MIDDLE(1,2);
    MZI_Equations(2,2) = EQ_MZI2_MIDDLE(2,2);
    MZI_Equations(3,2) = EQ_MZI2_MIDDLE(3,2); 
    EQ_MZI4_MIDDLE = Equation_Setting_Matrix_MZI2_After*EQ_DMZI4_H44;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    MZI_Equations(1,3) = EQ_MZI4_MIDDLE(1,3);
    MZI_Equations(2,3) = EQ_MZI4_MIDDLE(2,3);
    MZI_Equations(3,3) = EQ_MZI4_MIDDLE(3,3);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Type2 MZI3 

    Equation13 = MZI_Equations(3,1)*((0.5*X3)-0.5)+MZI_Equations(3,2)*((-0.5j*X3)-0.5j) == 0;

    MZI_Equations_Update = MZI_Equations;
    EQ_MZI3 = [1 -0.5j*X3-0.5j;1 0.5-(0.5*X3)];
    EQ_MZI3_H33 = [EQ_MZI3(1,1) EQ_MZI3(1,2) 0 ;EQ_MZI3(2,1) EQ_MZI3(2,2) 0 ;0 0 1];
    MZI_Equations_Update_Middle = MZI_Equations_Update*EQ_MZI3_H33;
    MZI_Equations_Update(1,2) = MZI_Equations_Update_Middle(1,2);
    MZI_Equations_Update(2,2) = MZI_Equations_Update_Middle(2,2);
    MZI_Equations_Update(3,2) = MZI_Equations_Update_Middle(3,2);     

    Equation15 =  MZI_Equations_Update(3,2)*((0.5*X5)-0.5)+MZI_Equations_Update(3,3)*((-0.5j*X5)-0.5j) == 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Type2 MZI5
    EQ_MZI5 = [1 -0.5j*X5-0.5j;1 0.5-(0.5*X5)];
    EQ_MZI5_H33 = [1 0 0;0 EQ_MZI5(1,1) EQ_MZI5(1,2);0 EQ_MZI5(2,1) EQ_MZI5(2,2)];
    MZI_Equations_Update_Middle = MZI_Equations_Update*EQ_MZI5_H33;    
    MZI_Equations_Update(1,3) = MZI_Equations_Update_Middle(1,3);
    MZI_Equations_Update(2,3) = MZI_Equations_Update_Middle(2,3);
    MZI_Equations_Update(3,3) = MZI_Equations_Update_Middle(3,3);  

    Equation16 = MZI_Equations_Update(3,3) - 1 == 0;
    Equation17 = MZI_Equations_Update(2,3) == 0;
    Equation18 = MZI_Equations_Update(1,3) == 0;
    
    [x1,x2,x3,x4,x5] = solve(Equation13,Equation15,Equation16,Equation17,Equation18,X1,X2,X3,X4,X5);

    x1 =  double(x1);
    x2 =  double(x2);
    x3 =  double(x3);
    x4 =  double(x4);
    x5 =  double(x5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    X1
    P_phi1 = 0:0.001:2*pi;
    N_Phi1 = length(P_phi1);
    N_P1 = length(x1);
    Judgement_Potential_Phi1 = zeros(N_P1,N_Phi1);
    Potential_Phi1_Min = zeros(1,N_P1);
    Potential_Phi1 = zeros(1,N_P1);
    for k = 1:N_P1
    for l = 1:N_Phi1
        Judgement_Potential_Phi1(k,l) = abs(exp(-1j*P_phi1(l))-x1(k));
    end
    Potential_Phi1_Min(k) = min(Judgement_Potential_Phi1(k,:));
    Position(k) = find(Judgement_Potential_Phi1(k,:) == Potential_Phi1_Min(k));
    Potential_Phi1(k) = P_phi1(Position(k));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    X2
    P_phi2 = 0:0.001:2*pi;
    N_Phi2 = length(P_phi2);
    N_P2 = length(x2);
    Judgement_Potential_Phi2 = zeros(N_P2,N_Phi2);
    Potential_Phi2_Min = zeros(1,N_P2);
    Potential_Phi2 = zeros(1,N_P2);    
    for k = 1:N_P2
    for l = 1:N_Phi2
        Judgement_Potential_Phi2(k,l) = abs(exp(-1j*P_phi2(l))-x2(k));
    end
    Potential_Phi2_Min(k) = min(Judgement_Potential_Phi2(k,:));
    Position(k) = find(Judgement_Potential_Phi2(k,:) == Potential_Phi2_Min(k));
    Potential_Phi2(k) = P_phi2(Position(k));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    X4
    P_phi4 = 0:0.001:2*pi;
    N_Phi4 = length(P_phi4);
    N_P4 = length(x4);
    Judgement_Potential_Phi4 = zeros(N_P4,N_Phi4);
    Potential_Phi4_Min = zeros(1,N_P4);
    Potential_Phi4 = zeros(1,N_P4);    
    for k = 1:N_P4
    for l = 1:N_Phi4
        Judgement_Potential_Phi4(k,l) = abs(exp(-1j*P_phi4(l))-x4(k));
    end
    Potential_Phi4_Min(k) = min(Judgement_Potential_Phi4(k,:));
    Position(k) = find(Judgement_Potential_Phi4(k,:) == Potential_Phi4_Min(k));
    Potential_Phi4(k) = P_phi4(Position(k));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    X3
    P_theta3 = 0:0.001:2*pi;
    N_theta3 = length(P_theta3);
    N_T3 = length(x3);
    Judgement_Potential_theta3 = zeros(N_T3,N_theta3);
    Potential_Theta3_Min = zeros(1,N_T3);
    Potential_Theta3 = zeros(1,N_T3);    
    for k = 1:N_T3
    for l = 1:N_theta3
        Judgement_Potential_theta3(k,l) = abs(exp(-1j*P_theta3(l))-x3(k));
    end
    Potential_Theta3_Min(k) = min(Judgement_Potential_theta3(k,:));
    Position(k) = find(Judgement_Potential_theta3(k,:) == Potential_Theta3_Min(k));
    Potential_Theta3(k) = P_theta3(Position(k));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    X5
    P_theta5 = 0:0.001:2*pi;
    N_theta5 = length(P_theta5);
    N_T5 = length(x5);
    Judgement_Potential_theta5 = zeros(N_T5,N_theta5);
    Potential_Theta5_Min = zeros(1,N_T5);
    Potential_Theta5 = zeros(1,N_T5);     
    for k = 1:N_T5
    for l = 1:N_theta5
        Judgement_Potential_theta5(k,l) = abs(exp(-1j*P_theta5(l))-x5(k));
    end
    Potential_Theta5_Min(k) = min(Judgement_Potential_theta5(k,:));
    Position(k) = find(Judgement_Potential_theta5(k,:) == Potential_Theta5_Min(k));
    Potential_Theta5(k) = P_theta5(Position(k));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------Rebuild------------------------
    MZI_Phi_1 = Potential_Phi1(Building_Parameter_Use_Num);
    MZI_Phi_2 = Potential_Phi2(Building_Parameter_Use_Num);
    MZI_Phi_4 = Potential_Phi4(Building_Parameter_Use_Num);
    MZI_Theta_3 = Potential_Theta3(Building_Parameter_Use_Num);
    MZI_Theta_5 = Potential_Theta5(Building_Parameter_Use_Num);

    [DMZI_Set1] = DMZI_Setup(MZI_Theta_1,MZI_Phi_1);
    [DMZI_Set2] = DMZI_Setup(MZI_Theta_2,MZI_Phi_2);
    [DMZI_Set4] = DMZI_Setup(MZI_Theta_4,MZI_Phi_4);

    DMZI_Set1_H44 = [DMZI_Set1(1,1) DMZI_Set1(1,2) 0 0;...
                     DMZI_Set1(2,1) DMZI_Set1(2,2) 0 0;...
                     0                 0           1 0;...
                     0                 0           0 1];
    DMZI_Set1_H44_Dagger = transpose(conj(DMZI_Set1_H44));  
    DMZI_Set2_H44 = [1                 0           0 0;...
                     0 DMZI_Set2(1,1) DMZI_Set2(1,2) 0;...
                     0 DMZI_Set2(2,1) DMZI_Set2(2,2) 0;...
                     0                 0           0 1];
    DMZI_Set2_H44_Dagger = transpose(conj(DMZI_Set2_H44));
    DMZI_Set4_H44 = [1 0           0                 0;...
                     0 1           0                 0;...
                     0 0 DMZI_Set4(1,1) DMZI_Set4(1,2);...
                     0 0 DMZI_Set4(2,1) DMZI_Set4(2,2)];
    DMZI_Set4_H44_Dagger = transpose(conj(DMZI_Set4_H44));

    Rebuild_Matrix = Venify_APP_TSU4*DMZI_Set1_H44_Dagger*DMZI_Set2_H44_Dagger*DMZI_Set4_H44_Dagger;

    DMZI_Set3_H44_Dagger = [0.5*exp(-1j*MZI_Theta_3)-0.5 -0.5j*exp(-1j*MZI_Theta_3)-0.5j 0 0;...
                            -0.5j*exp(-1j*MZI_Theta_3)-0.5j 0.5-0.5*exp(-1j*MZI_Theta_3) 0 0;...
                            0                                 0                          1 0;...
                            0                                 0                          0 1];
    DMZI_Set5_H44_Dagger = [1                                 0                          0 0;...
                            0 0.5*exp(-1j*MZI_Theta_5)-0.5 -0.5j*exp(-1j*MZI_Theta_5)-0.5j 0;...
                            0 -0.5j*exp(-1j*MZI_Theta_5)-0.5j 0.5-0.5*exp(-1j*MZI_Theta_5) 0;...
                            0                                 0                          0 1];
    %---------------Third procedure Equations Set up-------------------
%     syms P3 P5 T6 P6
    P3 = sym('P3');
    P5 = sym('P5');
    T6 = sym('T6');
    P6 = sym('P6');
    D_MZI_Dagger3 = [P3*DMZI_Set3_H44_Dagger(1,1) DMZI_Set3_H44_Dagger(1,2);...
                  P3*DMZI_Set3_H44_Dagger(2,1) DMZI_Set3_H44_Dagger(2,2)];
    D_MZI_Dagger5 = [P5*DMZI_Set5_H44_Dagger(2,2) DMZI_Set5_H44_Dagger(2,3);...
                  P5*DMZI_Set5_H44_Dagger(3,2) DMZI_Set5_H44_Dagger(3,3)];
    D_MZI_Dagger6 = [P6*(0.5*T6-0.5) -0.5j*T6-0.5j;...
                  P6*(-0.5j*T6-0.5j) 0.5-0.5*T6];
    D_MZI_Dagger3_H44 = [D_MZI_Dagger3(1,1) D_MZI_Dagger3(1,2) 0 0;...
                         D_MZI_Dagger3(2,1) D_MZI_Dagger3(2,2) 0 0;...
                         0                  0                  1 0;...
                         0                  0                  0 1];
    D_MZI_Dagger5_H44 = [1 0                  0                  0;...
                         0 D_MZI_Dagger5(1,1) D_MZI_Dagger5(1,2) 0;...
                         0 D_MZI_Dagger5(2,1) D_MZI_Dagger5(2,2) 0;...
                         0 0                  0                  1];
    D_MZI_Dagger6_H44 = [D_MZI_Dagger6(1,1) D_MZI_Dagger6(1,2) 0 0;...
                         D_MZI_Dagger6(2,1) D_MZI_Dagger6(2,2) 0 0;...
                         0                  0                  1 0;...
                         0                  0                  0 1];   
    IIIrd_Equations_Setup_A = Rebuild_Matrix*D_MZI_Dagger3_H44*D_MZI_Dagger5_H44;
    Equation20 = IIIrd_Equations_Setup_A(2,1)*(0.5*T6-0.5)+IIIrd_Equations_Setup_A(2,2)*(-0.5j*T6-0.5j) == 0;
    IIIrd_Equations_Setup_B = IIIrd_Equations_Setup_A*D_MZI_Dagger6_H44;
    Equation21 = IIIrd_Equations_Setup_B(2,2) == 1;
    Equation22 = IIIrd_Equations_Setup_B(1,2) == 0;
    Equation23 = IIIrd_Equations_Setup_B(1,1) == 1;
    [p3,p5,t6,p6] = solve(Equation20,Equation21,Equation22,Equation23,P3,P5,T6,P6); 
    p3 =  double(p3);
    p5 =  double(p5);
    t6 =  double(t6);
    p6 =  double(p6);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    P3
    P_phi3 = 0:0.001:2*pi;
    N_Phi3 = length(P_phi3);
    N_P3 = length(p3);
    Judgement_Potential_Phi3 = zeros(N_P3,N_Phi3);
    Potential_Phi3_Min = zeros(1,N_P3);
    Potential_Phi3 = zeros(1,N_P3);
    for k = 1:N_P3
    for l = 1:N_Phi3
        Judgement_Potential_Phi3(k,l) = abs(exp(-1j*P_phi3(l))-p3(k));
    end
    Potential_Phi3_Min(k) = min(Judgement_Potential_Phi3(k,:));
    Position(k) = find(Judgement_Potential_Phi3(k,:) == Potential_Phi3_Min(k));
    Potential_Phi3(k) = P_phi3(Position(k));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    P5
    P_phi5 = 0:0.001:2*pi;
    N_Phi5 = length(P_phi5);
    N_P5 = length(p5);
    Judgement_Potential_Phi5 = zeros(N_P5,N_Phi5);
    Potential_Phi5_Min = zeros(1,N_P5);
    Potential_Phi5 = zeros(1,N_P5);
    for k = 1:N_P5
    for l = 1:N_Phi5
        Judgement_Potential_Phi5(k,l) = abs(exp(-1j*P_phi5(l))-p5(k));
    end
    Potential_Phi5_Min(k) = min(Judgement_Potential_Phi5(k,:));
    Position(k) = find(Judgement_Potential_Phi5(k,:) == Potential_Phi5_Min(k));
    Potential_Phi5(k) = P_phi5(Position(k));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    P6
    P_phi6 = 0:0.001:2*pi;
    N_Phi6 = length(P_phi6);
    N_P6 = length(p6);
    Judgement_Potential_Phi6 = zeros(N_P6,N_Phi6);
    Potential_Phi6_Min = zeros(1,N_P6);
    Potential_Phi6 = zeros(1,N_P6);
    for k = 1:N_P6
    for l = 1:N_Phi6
        Judgement_Potential_Phi6(k,l) = abs(exp(-1j*P_phi6(l))-p6(k));
    end
    Potential_Phi6_Min(k) = min(Judgement_Potential_Phi6(k,:));
    Position(k) = find(Judgement_Potential_Phi6(k,:) == Potential_Phi6_Min(k));
    Potential_Phi6(k) = P_phi6(Position(k));
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------Result Transfer-----------------
    %                    X3
    P_theta6 = 0:0.001:2*pi;
    N_theta6 = length(P_theta6);
    N_T6 = length(t6);
    Judgement_Potential_theta6 = zeros(N_T6,N_theta6);
    Potential_Theta6_Min = zeros(1,N_T6);
    Potential_Theta6 = zeros(1,N_T6);    
    for k = 1:N_T6
    for l = 1:N_theta6
        Judgement_Potential_theta6(k,l) = abs(exp(-1j*P_theta6(l))-t6(k));
    end
    Potential_Theta6_Min(k) = min(Judgement_Potential_theta6(k,:));
    Position(k) = find(Judgement_Potential_theta6(k,:) == Potential_Theta6_Min(k));
    Potential_Theta6(k) = P_theta6(Position(k));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------Rebuild------------------------

    MZI_Phi_3 = Potential_Phi3(Building_Parameter_Use_Num_2);
    MZI_Phi_5 = Potential_Phi5(Building_Parameter_Use_Num_2);
    MZI_Phi_6 = Potential_Phi6(Building_Parameter_Use_Num_2);
    MZI_Theta_6 = Potential_Theta6(Building_Parameter_Use_Num_2);
    [DMZI_Set3] = DMZI_Setup(MZI_Theta_3,MZI_Phi_3);
    [DMZI_Set5] = DMZI_Setup(MZI_Theta_5,MZI_Phi_5);
    [DMZI_Set6] = DMZI_Setup(MZI_Theta_6,MZI_Phi_6); 

    DMZI_Set3_H44 = [DMZI_Set3(1,1) DMZI_Set3(1,2) 0 0;...
                     DMZI_Set3(2,1) DMZI_Set3(2,2) 0 0;...
                     0              0              1 0;...
                     0              0              0 1];
    DMZI_Set3_H44_Dagger = transpose(conj(DMZI_Set3_H44));                  
    DMZI_Set5_H44 = [1 0              0              0;...
                     0 DMZI_Set5(1,1) DMZI_Set5(1,2) 0;...
                     0 DMZI_Set5(2,1) DMZI_Set5(2,2) 0;...
                     0 0              0              1];
    DMZI_Set5_H44_Dagger = transpose(conj(DMZI_Set5_H44));                  
    DMZI_Set6_H44 = [DMZI_Set6(1,1) DMZI_Set6(1,2) 0 0;...
                     DMZI_Set6(2,1) DMZI_Set6(2,2) 0 0;...
                     0              0              1 0;...
                     0              0              0 1];
    DMZI_Set6_H44_Dagger = transpose(conj(DMZI_Set6_H44));
    Matrix_Final = Venify_APP_TSU4*DMZI_Set1_H44_Dagger*DMZI_Set2_H44_Dagger*DMZI_Set4_H44_Dagger*...
                                DMZI_Set3_H44_Dagger*DMZI_Set5_H44_Dagger*DMZI_Set6_H44_Dagger;
    Phi(1) = MZI_Phi_1;
    Phi(2) = MZI_Phi_2;
    Phi(3) = MZI_Phi_3;
    Phi(4) = MZI_Phi_4;
    Phi(5) = MZI_Phi_5;
    Phi(6) = MZI_Phi_6;
    Theta(1) = MZI_Theta_1;
    Theta(2) = MZI_Theta_2;
    Theta(3) = MZI_Theta_3;
    Theta(4) = MZI_Theta_4;
    Theta(5) = MZI_Theta_5;
    Theta(6) = MZI_Theta_6;                            
    disp('The Application SU(4) multiplex by MZIs result is');
    Venify_Final_result = det(round(abs(Matrix_Final)));
    disp(round(abs(Matrix_Final)));
    if(Venify_Final_result == 1)
        disp('Setting up Phis and Thetas are Effective');
        disp('     ------- Theta 1-6 -------');
        disp(Theta);
        disp('     ------- Phi 1-6 -------');
        disp(Phi);
        disp('   -------Calculation Finish-------');
    else
        disp('Setting up Phis and Thetas are not accurate')
        disp('   -------Calculation Finish-------');
    end  
end
    function [DMZI_Set] = DMZI_Setup(Theta_Set,Phi_Set)
    DMZI_Set = 1j*exp(1j*Theta_Set/2).*[exp(1j*Phi_Set)*sin(Theta_Set/2) exp(1j*Phi_Set)*cos(Theta_Set/2);...
    cos(Theta_Set/2) -1*sin(Theta_Set/2)];
    end
end