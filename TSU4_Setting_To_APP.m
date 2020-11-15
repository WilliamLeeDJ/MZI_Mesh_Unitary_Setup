function [TSU4] = TSU4_Setting_To_APP(Phi,Theta)

    DMZI1 = DMZI_Setup(Theta(1),Phi(1));
    DMZI2 = DMZI_Setup(Theta(2),Phi(2));
    DMZI3 = DMZI_Setup(Theta(3),Phi(3));
    DMZI4 = DMZI_Setup(Theta(4),Phi(4));
    DMZI5 = DMZI_Setup(Theta(5),Phi(5));
    DMZI6 = DMZI_Setup(Theta(6),Phi(6));
    
    DMZI1_H44 = DMZI_H44_Setup(DMZI1,1);
    DMZI2_H44 = DMZI_H44_Setup(DMZI2,2);
    DMZI3_H44 = DMZI_H44_Setup(DMZI3,1);
    DMZI4_H44 = DMZI_H44_Setup(DMZI4,3);
    DMZI5_H44 = DMZI_H44_Setup(DMZI5,2);
    DMZI6_H44 = DMZI_H44_Setup(DMZI6,1);
    
    TSU4 = DMZI6_H44*DMZI5_H44*DMZI4_H44*DMZI3_H44*DMZI2_H44*DMZI1_H44;
    
    disp('   ----- Attention! -----   ');
    disp('The Setup Transfer Matrix is');
    disp(TSU4);
    disp('   ----- Finish Cal! -----   ');
    function [DMZI_Set] = DMZI_Setup(Theta_Set,Phi_Set)
    DMZI_Set = 1j*exp(1j*Theta_Set/2).*[exp(1j*Phi_Set)*sin(Theta_Set/2) exp(1j*Phi_Set)*cos(Theta_Set/2);...
    cos(Theta_Set/2) -1*sin(Theta_Set/2)];
    end
    function [DMZI_H44] = DMZI_H44_Setup(DMZI_Set_H22,MZI_Top_Channel)
        if MZI_Top_Channel == 1
            DMZI_H44 = [...
                DMZI_Set_H22(1,1) DMZI_Set_H22(1,2) 0 0;...
                DMZI_Set_H22(2,1) DMZI_Set_H22(2,2) 0 0;...
                0                 0                 1 0;...
                0                 0                 0 1];
        end
        if MZI_Top_Channel == 2
             DMZI_H44 = [...
                1 0                 0                 0;...
                0 DMZI_Set_H22(1,1) DMZI_Set_H22(1,2) 0;...
                0 DMZI_Set_H22(2,1) DMZI_Set_H22(2,2) 0;...
                0                 0                 0 1];
        end
        if MZI_Top_Channel == 3
             DMZI_H44 = [...
                1 0 0                 0;...
                0 1 0                 0;...
                0 0 DMZI_Set_H22(1,1) DMZI_Set_H22(1,2);...
                0 0 DMZI_Set_H22(2,1) DMZI_Set_H22(2,2)];
        end       
    end
end