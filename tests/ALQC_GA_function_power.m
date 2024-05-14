function cost =ALQC_GA_function_power(A,G_temp,T_amb_temp,Q_int_day,Q_int_day_EL,Q_int_day_MP,m1_LW,m1_EL,m1_MP,flag_B1,flag_B2,flag_B3,Num_B1,Num_B2,Num_B3,T_g,T_g_EL,T_g_MP,ini_zone_B1,ini_zone_B2,ini_zone_B3,ini_tank_B1,ini_tank_B2,ini_tank_B3)
    cost=diff_yd(0.32)+diff_yd(0.35)+diff_yd(0.38); 
    function total=diff_yd(yd)
    %Theta=[0.000119899318008265;-0.0111979189569418;0.0124138222190044;0.00969803592160545;-0.0116749761637758;0.00200303523686252;-0.0141980836690010;0.0136211426786303;-0.000228564530296183;-0.899356144517810;0.0378607977479390];
    b1=5.237981727599351E-4;b2=0.12362070796193345;b3=-0.0841844009399956;
    c1=0.055357759742934484;c2=-0.034120519224108234;c3=0.004836274437410222;
    d1=0.003947362053623362;d2=0.03171560897057907;
    f1=-2.337002912084548E-4;
    a1=-0.25269378439996903;
    g=0.1755537160875479;
    a1_d=1-a1;a2_d=a1;

    Day=4;
    hour=0:1:Day*24;%hour
    Delta_t=1;%1s
    h=0:Delta_t/3600:Day*24;
    G_solar=interp1(hour,G_temp,h,'spline');%solar radiation W/m^2
    T_amb=interp1(hour,T_amb_temp,h,'spline');%outside temperature every second ℃
    t_final=3600*24*Day;%total time:second
    dt=1;%Sample once per second
    Tc=300;%the sampling time of 5min 
    n_final=t_final/dt+1;
    N_final=t_final/Tc+1;
    C_air=10^3;%thermal capacity of air(J/kg ℃)
    Cw=4182; %specific heat of water
    T_DCW=15;%the tap water temperature 
    alpha_sol=0.3;alpha_inf=0.8;
    alpha1=1;%coefficient for heat conductivity
    
    %%  %%%%%%%%%%%%%%%%%%%%%%%building 1:LW building%%%%%%%%%%%%%%%%%%%%%%%
    U_exwall=0.8;%W/m^2 K
    U_roof=0.47;
    U_floor=0.35; 
    U_win=2.8;
    A_win_1=(2.7+2.8+0.7+3.8+5.8)*flag_B1.^2;%m^2 lower floor!!!
    A_win_2=(4.6+0.8+1.2+2.4+3.8+0.2+2.1)*flag_B1.^2;%m^2 upper floor!!
    R_win_1=1./(A_win_1*U_win);%thermal resistance of the external window (K/W)
    R_win_2=1./(A_win_2*U_win);%thermal resistance of the external window (K/W)
    g_win=0.8;%Total solar heat transmittance
    q_air=7.3*alpha_inf;%Airtightness m^3/h m^2
    rho_air=1.293; %density of air(kg/m^3)
    aa=8.5*flag_B1;%m Num_B1
    bb=10.6*flag_B1;%m
    H_room=2.6*flag_B1;%height of the room m
    A_roof=aa.*bb;%m^2
    A_ground=A_roof;%1*Num_B1
    V1=A_ground.*H_room;%1*Num_B1
    V2=A_ground.*H_room;%1*Num_B1
    M_air_1=rho_air*V1;%mass of the indoor air
    M_air_2=rho_air*V2;%mass of the indoor air
    A_exwall_1=(aa+bb)*2.*H_room-A_win_1;%area of external wall (m^2)
    A_exwall_2=(aa+bb)*2.*H_room-A_win_2;%area of external wall (m^2)
    m_inf_1=q_air*rho_air*(A_exwall_1+A_win_1)/3600;%air flow rate due to infiltration (kg/s)
    m_inf_2=q_air*rho_air*(A_exwall_2+A_win_2)/3600;%air flow rate due to infiltration (kg/s)
    h_in=50;%heat transfer coefficient of the inside surface of the wall (W/(m2 K))
    h_out=60;%heat transfer coefficient of the outside surface of the wall (W/(m2 K))
    M_exwall_1=180*A_exwall_1;%mass of the external wall(kg) %(0.013*2960+0.15*625+0.05*130+0.009*1620)*A_exwall_1  %153.31*A_exwall_1
    M_exwall_2=180*A_exwall_2;%mass of the external wall(kg)
    %13mm Gypsum 2960kg/m^3, 150 mm wooden frame 625kg/m^3, 50 mm mineral wool 130kg/m^3, 9mm wind shield board 1620kg/m^3(glass)
    C_exwall=16200*A_exwall_1./M_exwall_1;%thermal capacity of the external wall (J/(kg K))
                                         %%(0.013*2960*950+0.15*625*1180+0.05*130*837+0.009*1620*840)*A_exwall_1/M_exwall_1;   164870*A_exwall_1/M_exwall_1
    %13mm Gypsum 950J/kg K, 150 mm wooden frame 1180J/kg K, 50 mm mineral wool 837J/kg K, 9mm wind shield board 840J/kg K(glass)
    R_exwall_1=1./(A_exwall_1*U_exwall)-1./(A_exwall_1*h_in)-1./(A_exwall_1*h_out);%thermal resistance of the external wall (K/W) 
    R_exwall_2=1./(A_exwall_2*U_exwall)-1./(A_exwall_2*h_in)-1./(A_exwall_2*h_out);%thermal resistance of the external wall (K/W) 
    Q_solar_1=zeros(3600*24*Day+1,Num_B1);Q_solar_2=zeros(3600*24*Day+1,Num_B1);
    for i=1:Num_B1
       Q_solar_1(:,i)=alpha_sol*g_win*A_win_1(i)*G_solar';
       Q_solar_2(:,i)=alpha_sol*g_win*A_win_2(i)*G_solar';
    end
    
    Q_EH=4000;
    Mw=300*flag_B1.^3; %mass of water in the tank kg
    T_zone_1=ones(3600*24*Day+1,1)*ini_zone_B1;%air temperature of zone
    T_zone_2=ones(3600*24*Day+1,1)*ini_zone_B1;%air temperature of zone
    T_exwallin_1=T_zone_1-0.3;T_exwallin_2=T_zone_2-0.3;
    T_exwallout_1=-1.1*ones(3600*24*Day+1,Num_B1);T_exwallout_2=-1.1*ones(3600*24*Day+1,Num_B1);%outside surface temperature of the external wall (℃)
    Q_HVAC_1=zeros(3600*24*Day+1,Num_B1);Q_HVAC_2=zeros(3600*24*Day+1,Num_B1);
    T_tank=ones(3600*24*Day+1,1)*ini_tank_B1;%water temperature of tank
    Q_ER_1=8000*ones(3600*24*Day+1,Num_B1);
    Q_ER_2=8000*ones(3600*24*Day+1,Num_B1);
    S_ER_1=zeros(3600*24*Day+1,Num_B1);% at the start-up,half the electricity radiator is off
    S_ER_2=zeros(3600*24*Day+1,Num_B1);% at the start-up,half the electricity radiator is off
    S_EH=zeros(3600*24*Day+1,Num_B1);% at the start-up,half the electric heater is off
    
    %%  %%%%%%%%%%%%%%%%%%%%%%%building 2:EL building%%%%%%%%%%%%%%%%%%%%%%%
    H_room_EL=2.5*flag_B2;%height of the room m
    U_exwall_EL=0.17;%W/m^2 K
    U_roof_EL=0.28;
    U_floor_EL=0.79; 
    U_door_EL=1;
    U_win_EL=0.7;
    A_win_EL=(3+4.5+3+3+3)*flag_B2.^2;%m^2
    A_door_EL=3*flag_B2.^2;
    R_win_EL=1./(A_win_EL*U_win_EL);%thermal resistance of the external window (K/W)
    R_door_EL=1./(A_door_EL*U_door_EL);%thermal resistance of the external door (K/W)
    g_win_EL=0.5;%Total solar heat transmittance
    q_air_EL=0.4*alpha_inf;%Airtightness m^3/h m^2
    rho_air_EL=1.293; %density of air(kg/m^3)
    cc=15*flag_B2;dd=10*flag_B2;
    A_roof_EL=cc.*dd;%m^2
    A_ground_EL=A_roof_EL;
    V_EL=A_ground_EL.*H_room_EL;
    M_air_EL=rho_air_EL*V_EL;%mass of the indoor air
    A_exwall_EL=(cc+dd)*2.*H_room_EL-A_win_EL-A_door_EL;
    m_inf_EL=q_air_EL*rho_air_EL*(A_exwall_EL+A_win_EL+A_door_EL+A_roof_EL)/3600;%kg/s
    h_in_EL=50;%heat transfer coefficient of the inside surface of the wall (W/(m2 K))
    h_out_EL=60;%heat transfer coefficient of the outside surface of the wall (W/(m2 K))
    M_exwall_EL=(0.012*530+0.212*11+0.025*970+0.022*1.293+0.028*1250)*A_exwall_EL;%0.012 m plasterboard 530kg/m^3, 0.212 m mineral wool 11kg/m^3, 0.025 m wind shield board 970kg/m^3, 0.022 m ventilation gap 1.293kg/m^3(air), 0.028 m cladding 1250kg/m^3
    C_exwall_EL=(0.012*530*1090+0.212*11*837+0.025*970*1090+0.022*1.293*1000+0.028*1250*1000)*A_exwall_EL./M_exwall_EL; %0.012 m plasterboard 1090 J/(kg K), 0.212 m mineral wool 837J/(kg K), 0.025 m wind shield board 1090J/kg K, 0.022 m ventilation gap 10^3J/kg K(air), 0.028 m cladding 1000J/kg K
    R_exwall_EL=1./(A_exwall_EL*U_exwall_EL)-1./(A_exwall_EL*h_in_EL)-1./(A_exwall_EL*h_out_EL);
    Q_solar_EL=zeros(3600*24*Day+1,Num_B2);
    for i=1:Num_B2
       Q_solar_EL(:,i)=alpha_sol*g_win_EL*A_win_EL(i)*G_solar';
    end
    % tank
    Q_EH_EL=5000;
    Q_GSHP_EL=7500;
    Mw_EL=300*flag_B2.^3; %mass of water in the tank kg
    d_EL=0.59*flag_B2;
    r_EL=d_EL/2;%internal diameter of the tank
    Atank_EL=pi*r_EL.^2;%cross-section area of the storage tank
    kw_EL=0.65;%(W/m K)thermal conductivity of water 50℃:0.64  destratification conductivity:34.6
    Htank_EL=0.925*flag_B2;%height of each layer (2 node)
    rho_case_EL=7850;%kg/m^3
    Vcase_EL=2*pi*Htank_EL.*((r_EL+0.0006).^2-r_EL.^2);
    Mcase_EL=rho_case_EL*Vcase_EL;
    Ccase_EL=502.416;%J/kg K
    hcase_EL=1500*flag_B2;%kcase*nu1*Htank;
    
    %space heating: hydronic radiator
    M_HR_EL=6;%mass of water inside the radiator(kg)
    A_HR_EL=7;%area of the radiator(m^2)
    T_HR_EL=40;%outlet temperature of the radiator
    m_HR_EL=0.15;%0.04 kg/(m^2 s)
    COP_GSHP_EL=3;
    
    T_zone_EL=ones(3600*24*Day+1,1)*ini_zone_B2;
    T_exwallin_EL=T_zone_EL-0.3;
    T_exwallout_EL=-1.1*ones(3600*24*Day+1,Num_B2);
    Q_HVAC_EL=zeros(3600*24*Day+1,Num_B2);
    T_tank_1_EL=ones(3600*24*Day+1,1)*ini_tank_B2;%the bottom layer of the tank
    T_tank_2_EL=ones(3600*24*Day+1,1)*ini_tank_B2;%the top layer of the tank
    S_EH_EL=zeros(3600*24*Day+1,Num_B2);
    S_GSHP_EL=zeros(3600*24*Day+1,Num_B2);
    Tcase_EL=ones(3600*24*Day+1,1)*ini_tank_B2;
    Q_HR_EL=30000*ones(3600*24*Day+1,Num_B2);S_HR_EL=zeros(3600*24*Day+1,Num_B2);
    T_w_EL=52*ones(3600*24*Day+1,Num_B2);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%building 3: MP building%%%%%%%%%%%%%%%%%%%%%%%
    H_room_MP=2.6*flag_B3;%height of the room m
    U_exwall_MP=0.08;%W/m^2 K
    U_roof_MP=0.07;
    U_floor_MP=0.08; 
%     U_door_MP=0.5;
    U_win_MP=0.8;
    A_win_1_MP=(2.7+2.8+0.7+3.8+5.8)*flag_B3.^2;%m^2 lower floor
    A_win_2_MP=(4.6+0.8+1.2+2.4+3.8+0.2+2.1)*flag_B3.^2;%m^2 upper floor
    A_door_MP=2.58*flag_B3.^2;
    R_win_1_MP=1./(A_win_1_MP*U_win_MP);%thermal resistance of the external window (K/W)
    R_win_2_MP=1./(A_win_2_MP*U_win_MP);%thermal resistance of the external window (K/W)
%     R_door_MP=1./(A_door_MP*U_door_MP);%thermal resistance of the external door (K/W)
    g_win_MP=0.5;%Total solar heat transmittance
    q_air_MP=0.7*alpha_inf;%Airtightness m^3/h m^2
    rho_air_MP=1.293; %density of air(kg/m^3)
    ee=8.5*flag_B3;ff=10.6*flag_B3;
    A_roof_MP=ee.*ff;%m^2
    A_ground_MP=A_roof_MP;
    V_MP=A_ground_MP.*H_room_MP;
    M_air_MP=rho_air_MP*V_MP;%mass of the indoor air
    A_exwall_1_MP=(ee+ff)*2.*H_room_MP-A_win_1_MP-A_door_MP;
    A_exwall_2_MP=(ee+ff)*2.*H_room_MP-A_win_2_MP;
    m_inf_1_MP=q_air_MP*rho_air_MP*(A_exwall_1_MP+A_win_1_MP+A_door_MP)/3600;%kg/s
    m_inf_2_MP=q_air_MP*rho_air_MP*(A_exwall_2_MP+A_win_2_MP+A_roof_MP)/3600;%kg/s
    h_in_MP=6;%heat transfer coefficient of the inside surface of the wall (W/(m2 K))
    h_out_MP=8;%heat transfer coefficient of the outside surface of the wall (W/(m2 K))
    M_exwall_1_MP=(0.13*1500+0.34*30+0.09*1500)*A_exwall_1_MP;
    M_exwall_2_MP=(0.13*1500+0.34*30+0.09*1500)*A_exwall_2_MP;
    %0.13 m Light weight concrete block 1500kg/m^3, 0.34m  polyurethane 30kg/m^3, 0.09 m Light weight concrete block 1500kg/m^3
    C_exwall_MP=(0.13*1500*879+0.34*30*1500+0.09*1500*879)*A_exwall_1_MP./M_exwall_1_MP; 
    %0.13 m Light weight concrete block 879J/(kg K), 0.34m  polyurethane 1500J/(kg K), 0.09 m Light weight concrete block 879J/(kg K) 
    R_exwall_1_MP=1./(A_exwall_1_MP*U_exwall_MP)-1./(A_exwall_1_MP*h_in_MP)-1./(A_exwall_1_MP*h_out_MP);
    R_exwall_2_MP=1./(A_exwall_2_MP*U_exwall_MP)-1./(A_exwall_2_MP*h_in_MP)-1./(A_exwall_2_MP*h_out_MP);
    Q_solar_1_MP=zeros(3600*24*Day+1,Num_B3);Q_solar_2_MP=zeros(3600*24*Day+1,Num_B3);
    for i=1:Num_B3
       Q_solar_1_MP(:,i)=alpha_sol*g_win_MP*A_win_1_MP(i)*G_solar';
       Q_solar_2_MP(:,i)=alpha_sol*g_win_MP*A_win_2_MP(i)*G_solar';
    end
    % tank
    Q_EH_MP=6000;
    Q_GSHP_MP=8000;COP_GSHP_MP=4.5;
    Mw_MP=500*flag_B3.^3; %mass of water in the tank kg
    d_MP=0.59*flag_B3;r_MP=d_MP/2;%internal diameter of the tank
    Atank_MP=pi*r_MP.^2;%cross-section area of the storage tank
    kw_MP=0.65;%(W/m K)thermal conductivity of water 50℃:0.64  destratification conductivity:34.6
    Htank_MP=0.925*flag_B3;%hight of each layer (2 node)
    rho_case_MP=7850;%kg/m^3
    Vcase_MP=2*pi*Htank_MP.*((r_MP+0.0006).^2-r_MP.^2);
    Mcase_MP=rho_case_MP*Vcase_MP;
    Ccase_MP=502.416;%J/kg K
    hcase_MP=1500*flag_B3;%kcase*nu1*Htank;
    
    %space heating: hydronic radiator
    M_HR_1_MP=4.5;M_HR_2_MP=4.5;%mass of water inside the radiator(kg)
    A_HR_1_MP=5;A_HR_2_MP=5;%area of the radiator(m^2)
    T_HR_1_MP=40;T_HR_2_MP=40;%outlet temperature of the radiator
    m_HR_MP=0.07*ones(2,1);%0.04 kg/(m^2 s)
    
    
    T_zone_1_MP=ones(3600*24*Day+1,1)*ini_zone_B3;T_zone_2_MP=ones(3600*24*Day+1,1)*ini_zone_B3;
    T_exwallin_1_MP=T_zone_1_MP-0.3;T_exwallin_2_MP=T_zone_2_MP-0.3;
    T_exwallout_1_MP=-1.1*ones(3600*24*Day+1,Num_B3);T_exwallout_2_MP=-1.1*ones(3600*24*Day+1,Num_B3);
    Q_HVAC_1_MP=zeros(3600*24*Day+1,Num_B3);Q_HVAC_2_MP=zeros(3600*24*Day+1,Num_B3);
    T_tank_1_MP=ones(3600*24*Day+1,1)*ini_tank_B3;
    T_tank_2_MP=ones(3600*24*Day+1,1)*ini_tank_B3;
    S_EH_MP=zeros(3600*24*Day+1,Num_B3);
    S_GSHP_MP=zeros(3600*24*Day+1,Num_B3);
    Tcase_MP=ones(3600*24*Day+1,1)*ini_tank_B3;
    Q_HR_1_MP=15000*ones(3600*24*Day+1,Num_B3);Q_HR_2_MP=15000*ones(3600*24*Day+1,Num_B3);
    S_HR_1_MP=zeros(3600*24*Day+1,Num_B3);S_HR_2_MP=zeros(3600*24*Day+1,Num_B3);
    T_w_1_MP=45*ones(3600*24*Day+1,Num_B3);T_w_2_MP=45*ones(3600*24*Day+1,Num_B3);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%% simulation 1: yd=0.32  %%%%%%%%%%%%%%%%%%%%%%%
    U1=zeros(N_final,1);%the setpoint shift of indoor temperature
    U2=zeros(N_final,1);%the setpoint shift of tank temperature
    delta_U1=zeros(N_final,1);%delta_U1(n_step)=U1(n_step)-U1(n_step-1)
    delta_U2=zeros(N_final,1);%delta_U2(n_step)=U2(n_step)-U2(n_step-1)
    delta_G_solar=zeros(N_final,1);
    delta_T_amb=zeros(N_final,1);
    xd=zeros(9,N_final);
    %J=0;
    T_set=zeros(N_final,1);
    T_set_tank=zeros(N_final,1);
    
    p_total=zeros(3600*24*Day+1,1);%the aggregated average power
    p_norm=zeros(3600*24*Day+1,1);%the normalized aggregated average power
    p_total_LW=zeros(3600*24*Day+1,1);%the aggregated average power
    p_total_EL=zeros(3600*24*Day+1,1);%the aggregated average power
    p_total_MP=zeros(3600*24*Day+1,1);%the aggregated average power
    
    P_norm=zeros(N_final,1);%the normalized average power of last 5 min
    y=zeros(N_final,1);
    e=zeros(N_final,1);
    Q=diag([0 0 0 0 0 0 0 1 1]);
    R=diag([0.10342775960910355 0.0499945377004648]);
    P=Q;
    %%LS
    beta=10^A(1);alpha=10^A(2);sigma=10^A(3);
    z=zeros(1,N_final); %SPM: z=theta'*Phi
    theta1=b1*ones(N_final,1);theta2=b2*ones(N_final,1);theta3=b3*ones(N_final,1);
    theta4=c1*ones(N_final,1);theta5=c2*ones(N_final,1);theta6=c3*ones(N_final,1);
    theta7=d1*ones(N_final,1);theta8=d2*ones(N_final,1);
    theta9=f1*ones(N_final,1);
    theta10=a1_d*ones(N_final,1); theta11=a2_d*ones(N_final,1);
    epsilon=zeros(N_final,1);%estimation error
    % 定义对角线上的元素值
    D = [10^A(4), 10^A(5), 10^A(6), 10^A(6), 10^A(6), 10^A(7), 10^A(7), 10^A(6), 10^A(4), 10^A(8), 10^A(5)];
    % 创建对角矩阵P1
    P1 = diag(D);
    % 创建一个包含所有 theta 向量的 cell array
    thetas = {theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10,theta11};
    % 定义参数的上下界
    LB1 = [0.00005, 0.05, -0.12, 0.02, -0.08, 0.001, 0.001, 0.01, -0.005,1.1, -0.7];
    UB1 = [0.0009, 0.2, -0.04, 0.08, -0.01, 0.01, 0.09, 0.08, -0.0001,1.7, -0.1];

    for k=1:dt:n_final-1
        j=ceil(k/3600);
        %the number of internal heat gain
        int_num=mod(j,24);
        if int_num==0
           int_num=1;
        end
        n_step=ceil(k/Tc);
        T_set(n_step)=22+U1(n_step);
        T_set_tank(n_step)=64+U2(n_step);
    
    %%   %%%%%%%Building 1:LW building %%%%%%%
        for p1=1:Num_B1
            if T_amb(k+1)<8
               Q_ER_1(k+1,p1)=15000;Q_ER_2(k+1,p1)=15000;
            end
            %zone1:lower floor
            if T_zone_1(k,p1)>T_set(n_step)+0.5
                S_ER_1(k+1,p1)=0;
            elseif T_zone_1(k,p1)<=T_set(n_step)-0.5
                S_ER_1(k+1,p1)=1;
            else
                S_ER_1(k+1,p1)=S_ER_1(k,p1);
            end
            T_exwallin_1(k+1,p1)=T_exwallin_1(k,p1)+2*(h_in*A_exwall_1(p1)*(T_zone_1(k,p1)-T_exwallin_1(k,p1))+(T_exwallout_1(k,p1)-T_exwallin_1(k,p1))/R_exwall_1(p1))*Delta_t/(M_exwall_1(p1)*C_exwall(p1));        
            T_exwallout_1(k+1,p1)=T_exwallout_1(k,p1)+2*(h_out*A_exwall_1(p1)*(T_amb(k)-T_exwallout_1(k,p1))+(T_exwallin_1(k,p1)-T_exwallout_1(k,p1))/R_exwall_1(p1)+0.3*0.65*G_solar(k)*A_exwall_1(p1))*Delta_t/(M_exwall_1(p1)*C_exwall(p1));        
            T_zone_1(k+1,p1)=T_zone_1(k,p1)+(h_in*A_exwall_1(p1)*(T_exwallin_1(k,p1)-T_zone_1(k,p1))+(T_amb(k)-T_zone_1(k,p1))/R_win_1(p1)+m_inf_1(p1)*C_air*(T_amb(k)-T_zone_1(k,p1))+A_ground(p1)*U_floor*(T_g(p1)-T_zone_1(k,p1))+A_ground(p1)*U_floor*(T_zone_2(k,p1)-T_zone_1(k,p1))+Q_HVAC_1(k,p1)+Q_int_day(int_num,p1)/2+Q_solar_1(k,p1))*Delta_t/(M_air_1(p1)*C_air);        
            Q_HVAC_1(k+1,p1)=Q_ER_1(k+1,p1)*S_ER_1(k+1,p1);
            %zone2:upper floor
            if T_zone_2(k,p1)>T_set(n_step)+0.5
                S_ER_2(k+1,p1)=0;
            elseif T_zone_2(k,p1)<=T_set(n_step)-0.5
                S_ER_2(k+1,p1)=1;
            else
                S_ER_2(k+1,p1)=S_ER_2(k,p1);
            end
            T_exwallin_2(k+1,p1)=T_exwallin_2(k,p1)+2*(h_in*A_exwall_2(p1)*(T_zone_2(k,p1)-T_exwallin_2(k,p1))+(T_exwallout_2(k,p1)-T_exwallin_2(k,p1))/R_exwall_2(p1))*Delta_t/(M_exwall_2(p1)*C_exwall(p1));
            T_exwallout_2(k+1,p1)=T_exwallout_2(k,p1)+2*(h_out*A_exwall_2(p1)*(T_amb(k)-T_exwallout_2(k,p1))+(T_exwallin_2(k,p1)-T_exwallout_2(k,p1))/R_exwall_2(p1)+0.3*0.65*G_solar(k)*A_exwall_2(p1))*Delta_t/(M_exwall_2(p1)*C_exwall(p1));
            T_zone_2(k+1,p1)=T_zone_2(k,p1)+(h_in*A_exwall_2(p1)*(T_exwallin_2(k,p1)-T_zone_2(k,p1))+(T_amb(k)-T_zone_2(k,p1))/R_win_2(p1)+m_inf_2(p1)*C_air*(T_amb(k)-T_zone_2(k,p1))+A_roof(p1)*U_roof*(T_amb(k)-T_zone_2(k,p1))+A_ground(p1)*U_floor*(T_zone_1(k,p1)-T_zone_2(k,p1))+Q_HVAC_2(k,p1)+Q_int_day(int_num,p1)/2+Q_solar_2(k,p1))*Delta_t/(M_air_2(p1)*C_air);
            Q_HVAC_2(k+1,p1)=Q_ER_2(k+1,p1)*S_ER_2(k+1,p1);
    
            %tank control
            if T_tank(k,p1)<T_set_tank(n_step)-1
                S_EH(k+1,p1)=1;
            elseif T_tank(k,p1)>T_set_tank(n_step)+1
                S_EH(k+1,p1)=0;
            else
                S_EH(k+1,p1)=S_EH(k,p1);
            end
            m1_prime_s2=m1_LW(p1,int_num)*(55-T_DCW)/(3600*(T_tank(k,p1)-T_DCW));
            Q_loss_s2=0.5*(T_tank(k,p1)-T_zone_1(k,p1));
            T_tank(k+1,p1)=T_tank(k,p1)+(S_EH(k,p1)*Q_EH+m1_prime_s2*Cw*(T_DCW-T_tank(k,p1))-Q_loss_s2)*Delta_t/(Mw(p1)*Cw);
            p_total_LW(k+1)=p_total_LW(k+1)+Q_ER_1(k+1,p1)*S_ER_1(k+1,p1)+Q_ER_2(k+1,p1)*S_ER_2(k+1,p1)+Q_EH*S_EH(k+1,p1);     
        end
    %%  %%%%%%%Building 2:EL building %%%%%%%
        for p2=1:Num_B2
            if T_zone_EL(k,p2)<=T_set(n_step)-0.5
                S_HR_EL(k+1,p2)=1;
            elseif T_zone_EL(k,p2)>T_set(n_step)+0.5
                S_HR_EL(k+1,p2)=0;
            else
                S_HR_EL(k+1,p2)=S_HR_EL(k,p2);
            end
            K_HR_EL=A_HR_EL*(T_w_EL(k,p2)-T_zone_EL(k,p2))^1.3;
            T_w_EL(k+1,p2)=T_w_EL(k,p2)+(m_HR_EL*S_HR_EL(k,p2)*Cw*(T_tank_1_EL(k,p2)-T_HR_EL)-Q_HR_EL(k,p2))*Delta_t/(M_HR_EL*Cw);
            T_exwallin_EL(k+1,p2)=T_exwallin_EL(k,p2)+2*(h_in_EL*A_exwall_EL(p2)*(T_zone_EL(k,p2)-T_exwallin_EL(k,p2))+(T_exwallout_EL(k,p2)-T_exwallin_EL(k,p2))/R_exwall_EL(p2))*Delta_t/(M_exwall_EL(p2)*C_exwall_EL(p2));
            T_exwallout_EL(k+1,p2)=T_exwallout_EL(k,p2)+2*(h_out_EL*A_exwall_EL(p2)*(T_amb(k)-T_exwallout_EL(k,p2))+(T_exwallin_EL(k,p2)-T_exwallout_EL(k,p2))/R_exwall_EL(p2)+0.3*0.65*G_solar(k)*A_exwall_EL(p2))*Delta_t/(M_exwall_EL(p2)*C_exwall_EL(p2));
            T_zone_EL(k+1,p2)=T_zone_EL(k,p2)+(h_in_EL*A_exwall_EL(p2)*(T_exwallin_EL(k,p2)-T_zone_EL(k,p2))+(T_amb(k)-T_zone_EL(k,p2))/R_win_EL(p2)+(T_amb(k)-T_zone_EL(k,p2))/R_door_EL(p2)+m_inf_EL(p2)*C_air*(T_amb(k)-T_zone_EL(k,p2))+A_roof_EL(p2)*U_roof_EL*(T_amb(k)-T_zone_EL(k,p2))+A_ground_EL(p2)*U_floor_EL*(T_g_EL(p2)-T_zone_EL(k,p2))+Q_int_day_EL(int_num,p2)+Q_HVAC_EL(k,p2)+Q_solar_EL(k,p2))*Delta_t/(M_air_EL(p2)*C_air);%
            Q_HR_EL(k+1,p2)=K_HR_EL*(T_w_EL(k+1,p2)-T_zone_EL(k+1,p2))*S_HR_EL(k+1,p2);
            Q_HVAC_EL(k+1,p2)=Q_HR_EL(k+1,p2);
    
            %tank control
            if (T_tank_1_EL(k,p2)<T_set_tank(n_step)-1)&&(T_tank_2_EL(k,p2)<90)
                S_GSHP_EL(k+1,p2)=1;
            elseif (T_tank_1_EL(k,p2)>T_set_tank(n_step)+1)||(T_tank_2_EL(k,p2)>=90)
                S_GSHP_EL(k+1,p2)=0;
            else
                S_GSHP_EL(k+1,p2)=S_GSHP_EL(k,p2);
            end
            if T_tank_1_EL(k,p2)<55 
                S_EH_EL(k+1,p2)=1;
            else
                S_EH_EL(k+1,p2)=0;
            end
            m1_prime_EL=m1_EL(p2,int_num)*(55-T_DCW)/(3600*(T_tank_1_EL(k,p2)-T_DCW));
            Q1_EL=alpha1*Atank_EL(p2)*kw_EL*(T_tank_2_EL(k,p2)-T_tank_1_EL(k,p2))/Htank_EL(p2);%convection between node 1 and 2
            Qc1_EL=hcase_EL(p2)*(pi*d_EL(p2)*Htank_EL(p2))*(Tcase_EL(k,p2)-T_tank_1_EL(k,p2));
            T_tank_1_EL(k+1,p2)=T_tank_1_EL(k,p2)+(S_EH_EL(k,p2)*Q_EH_EL+(m1_prime_EL+m_HR_EL*S_HR_EL(k,p2))*Cw*(T_tank_2_EL(k,p2)-T_tank_1_EL(k,p2))+Q1_EL+Qc1_EL)*Delta_t/(0.5*Mw_EL(p2)*Cw);
    
            Q2_EL=alpha1*Atank_EL(p2)*kw_EL*(T_tank_1_EL(k,p2)-T_tank_2_EL(k,p2))/Htank_EL(p2);%convection between node 1 and 2
            Qc2_EL=hcase_EL(p2)*(pi*d_EL(p2)*Htank_EL(p2))*(Tcase_EL(k,p2)-T_tank_2_EL(k,p2));
            T_tank_2_EL(k+1,p2)=T_tank_2_EL(k,p2)+(S_GSHP_EL(k,p2)*Q_GSHP_EL+m1_prime_EL*Cw*(T_DCW-T_tank_2_EL(k,p2))+m_HR_EL*S_HR_EL(k,p2)*Cw*(T_HR_EL-T_tank_2_EL(k,p2))+Q2_EL+Qc2_EL)*Delta_t/(0.5*Mw_EL(p2)*Cw);
    
            Qloss_EL=0.5*(Tcase_EL(k,p2)-T_zone_EL(k,p2));
            Tcase_EL(k+1,p2)=Tcase_EL(k,p2)-(Qc1_EL+Qc2_EL+Qloss_EL)*Delta_t/(Ccase_EL*Mcase_EL(p2));%case temperature of the tank (℃)
            p_total_EL(k+1)=p_total_EL(k+1)+Q_GSHP_EL*S_GSHP_EL(k+1,p2)/(COP_GSHP_EL)+Q_EH_EL*S_EH_EL(k+1,p2);    
        end
    %%   %%%%%%%Building 3:MP building %%%%%%% 
        for p3=1:Num_B3
            %zone1:lower floor
            if T_zone_1_MP(k,p3)<=T_set(n_step)-0.5
                S_HR_1_MP(k+1,p3)=1;
            elseif (T_zone_1_MP(k,p3)>T_set(n_step)+0.5)
                S_HR_1_MP(k+1,p3)=0;
            else
                S_HR_1_MP(k+1,p3)=S_HR_1_MP(k,p3);
            end
            K_HR_MP=A_HR_1_MP*(T_w_1_MP(k,p3)-T_zone_1_MP(k,p3))^1.3;
            T_w_1_MP(k+1,p3)=T_w_1_MP(k,p3)+(m_HR_MP(1)*S_HR_1_MP(k,p3)*Cw*(T_tank_1_MP(k,p3)-T_HR_1_MP)-Q_HR_1_MP(k,p3))*Delta_t/(M_HR_1_MP*Cw);
            T_exwallin_1_MP(k+1,p3)=T_exwallin_1_MP(k,p3)+2*(h_in_MP*A_exwall_1_MP(p3)*(T_zone_1_MP(k,p3)-T_exwallin_1_MP(k,p3))+(T_exwallout_1_MP(k,p3)-T_exwallin_1_MP(k,p3))/R_exwall_1_MP(p3))*Delta_t/(M_exwall_1_MP(p3)*C_exwall_MP(p3));
            T_exwallout_1_MP(k+1,p3)=T_exwallout_1_MP(k,p3)+2*(h_out_MP*A_exwall_1_MP(p3)*(T_amb(k)-T_exwallout_1_MP(k,p3))+(T_exwallin_1_MP(k,p3)-T_exwallout_1_MP(k,p3))/R_exwall_1_MP(p3)+0.3*0.65*G_solar(k)*A_exwall_1_MP(p3))*Delta_t/(M_exwall_1_MP(p3)*C_exwall_MP(p3));
            T_zone_1_MP(k+1,p3)=T_zone_1_MP(k,p3)+(h_in_MP*A_exwall_1_MP(p3)*(T_exwallin_1_MP(k,p3)-T_zone_1_MP(k,p3))+(T_amb(k)-T_zone_1_MP(k,p3))/R_win_1_MP(p3)+m_inf_1_MP(p3)*C_air*(T_amb(k)-T_zone_1_MP(k,p3))+A_ground_MP(p3)*U_floor_MP*(T_g_MP(p3)-T_zone_1_MP(k,p3))+A_ground_MP(p3)*U_floor_MP*(T_zone_2_MP(k,p3)-T_zone_1_MP(k,p3))+Q_int_day_MP(int_num,p3)/2+Q_solar_1_MP(k,p3)+Q_HVAC_1_MP(k,p3))*Delta_t/(M_air_MP(p3)*C_air);%
            Q_HR_1_MP(k+1,p3)=K_HR_MP*(T_w_1_MP(k+1,p3)-T_zone_1_MP(k+1,p3))*S_HR_1_MP(k+1,p3);
            Q_HVAC_1_MP(k+1,p3)=Q_HR_1_MP(k+1,p3);
    
            %zone2:upper floor
            if T_zone_2_MP(k,p3)<=T_set(n_step)-0.5
                S_HR_2_MP(k+1,p3)=1;
            elseif (T_zone_2_MP(k,p3)>T_set(n_step)+0.5)
                S_HR_2_MP(k+1,p3)=0;
            else
                S_HR_2_MP(k+1,p3)=S_HR_2_MP(k,p3);
            end
            K_HR_MP=A_HR_2_MP*(T_w_2_MP(k,p3)-T_zone_2_MP(k,p3))^1.3;
            T_w_2_MP(k+1,p3)=T_w_2_MP(k,p3)+(m_HR_MP(2)*S_HR_2_MP(k,p3)*Cw*(T_tank_2_MP(k,p3)-T_HR_2_MP)-Q_HR_2_MP(k,p3))*Delta_t/(M_HR_2_MP*Cw);
            T_exwallin_2_MP(k+1,p3)=T_exwallin_2_MP(k,p3)+2*(h_in_MP*A_exwall_2_MP(p3)*(T_zone_2_MP(k,p3)-T_exwallin_2_MP(k,p3))+(T_exwallout_2_MP(k,p3)-T_exwallin_2_MP(k,p3))/R_exwall_2_MP(p3))*Delta_t/(M_exwall_2_MP(p3)*C_exwall_MP(p3));
            T_exwallout_2_MP(k+1,p3)=T_exwallout_2_MP(k,p3)+2*(h_out_MP*A_exwall_2_MP(p3)*(T_amb(k)-T_exwallout_2_MP(k,p3))+(T_exwallin_2_MP(k,p3)-T_exwallout_2_MP(k,p3))/R_exwall_2_MP(p3)+0.3*0.65*G_solar(k)*A_exwall_2_MP(p3))*Delta_t/(M_exwall_2_MP(p3)*C_exwall_MP(p3));
            T_zone_2_MP(k+1,p3)=T_zone_2_MP(k,p3)+(h_in_MP*A_exwall_2_MP(p3)*(T_exwallin_2_MP(k,p3)-T_zone_2_MP(k,p3))+(T_amb(k)-T_zone_2_MP(k,p3))/R_win_2_MP(p3)+m_inf_2_MP(p3)*C_air*(T_amb(k)-T_zone_2_MP(k,p3))+A_roof_MP(p3)*U_roof_MP*(T_amb(k)-T_zone_2_MP(k,p3))+A_ground_MP(p3)*U_floor_MP*(T_zone_1_MP(k,p3)-T_zone_2_MP(k,p3))+Q_int_day_MP(int_num,p3)/2+Q_solar_2_MP(k,p3)+Q_HVAC_2_MP(k,p3))*Delta_t/(M_air_MP(p3)*C_air);%
            Q_HR_2_MP(k+1,p3)=K_HR_MP*(T_w_2_MP(k+1,p3)-T_zone_2_MP(k+1,p3))*S_HR_2_MP(k+1,p3);
            Q_HVAC_2_MP(k+1,p3)=Q_HR_2_MP(k+1,p3);
           %tank control
            if (T_tank_1_MP(k,p3)<T_set_tank(n_step)-1)&&(T_tank_2_MP(k,p3)<90)
                S_GSHP_MP(k+1,p3)=1;
            elseif (T_tank_1_MP(k,p3)>T_set_tank(n_step)+1)||(T_tank_2_MP(k,p3)>=90)
                S_GSHP_MP(k+1,p3)=0;
            else
                S_GSHP_MP(k+1,p3)=S_GSHP_MP(k,p3);
            end
            if T_tank_1_MP(k,p3)<55 
                S_EH_MP(k+1,p3)=1;
            else
                S_EH_MP(k+1,p3)=0;
            end
            m1_prime_MP=m1_MP(p3,int_num)*(55-T_DCW)/(3600*(T_tank_1_MP(k,p3)-T_DCW));
            Q1_MP=alpha1*Atank_MP(p3)*kw_MP*(T_tank_2_MP(k,p3)-T_tank_1_MP(k,p3))/Htank_MP(p3);%convection between node 1 and 2
            Qc1_MP=hcase_MP(p3)*(pi*d_MP(p3)*Htank_MP(p3))*(Tcase_MP(k,p3)-T_tank_1_MP(k,p3));
            T_tank_1_MP(k+1,p3)=T_tank_1_MP(k,p3)+(S_EH_MP(k,p3)*Q_EH_MP+(m1_prime_MP+m_HR_MP(1)*S_HR_1_MP(k,p3)+m_HR_MP(2)*S_HR_2_MP(k,p3))*Cw*(T_tank_2_MP(k,p3)-T_tank_1_MP(k,p3))+Q1_MP+Qc1_MP)*Delta_t/(0.5*Mw_MP(p3)*Cw);
    
            Q2_MP=alpha1*Atank_MP(p3)*kw_MP*(T_tank_1_MP(k,p3)-T_tank_2_MP(k,p3))/Htank_MP(p3);%convection between node 1 and 2
            Qc2_MP=hcase_MP(p3)*(pi*d_MP(p3)*Htank_MP(p3))*(Tcase_MP(k,p3)-T_tank_2_MP(k,p3));
            T_tank_2_MP(k+1,p3)=T_tank_2_MP(k,p3)+(S_GSHP_MP(k,p3)*Q_GSHP_MP+m1_prime_MP*Cw*(T_DCW-T_tank_2_MP(k,p3))+m_HR_MP(1)*S_HR_1_MP(k,p3)*Cw*T_HR_1_MP+m_HR_MP(2)*S_HR_2_MP(k,p3)*Cw*T_HR_2_MP-(m_HR_MP(1)*S_HR_1_MP(k,p3)+m_HR_MP(2)*S_HR_2_MP(k,p3))*Cw*T_tank_2_MP(k,p3)+Q2_MP+Qc2_MP)*Delta_t/(0.5*Mw_MP(p3)*Cw);
    
            Qloss_MP=0.5*(Tcase_MP(k,p3)-T_zone_1_MP(k,p3));
            Tcase_MP(k+1,p3)=Tcase_MP(k,p3)-(Qc1_MP+Qc2_MP+Qloss_MP)*Delta_t/(Ccase_MP*Mcase_MP(p3));%case temperature of the tank (℃)
            p_total_MP(k+1)=p_total_MP(k+1)+Q_GSHP_MP*S_GSHP_MP(k+1,p3)/(COP_GSHP_MP)+Q_EH_MP*S_EH_MP(k+1,p3);    
        end

        p_total(k+1)=p_total_LW(k+1)+p_total_EL(k+1)+p_total_MP(k+1);
        p_norm(k+1)=p_total(k+1)/((15+15+4)*1000*Num_B1+(5+7.5/3)*1000*Num_B2+(6+8/4.5)*1000*Num_B3);%GSHP
%% control law      
        if (mod(k-1,Tc)==0) &&(k<n_final-1)&&(k>1)
            P_norm(n_step)=sum(p_norm(k-(Tc-1):k))/Tc; %the normalized average power of 5 min
            y(n_step)=P_norm(n_step);    
            e(n_step)=y(n_step)-yd;

            %parameter LS algorithm
            z(n_step)=e(n_step);
            delta_G_solar(n_step)=G_solar(n_step)-G_solar(n_step-1);
            delta_T_amb(n_step)=T_amb(n_step)-T_amb(n_step-1);
            Phi=[delta_U1(max(n_step-1,1)),delta_U1(max(n_step-2,1)),delta_U1(max(n_step-3,1)),delta_U2(max(n_step-1,1)),delta_U2(max(n_step-2,1)),delta_U2(max(n_step-3,1)),delta_T_amb(max(n_step-1,1)),delta_T_amb(max(n_step-2,1)),delta_G_solar(max(n_step-1,1)),e(max(n_step-1,1)),e(max(n_step-2,1))]';    
            P1=1/beta.*(P1-P1*(Phi*Phi')*P1/((1+alpha*(Phi'*Phi)+(Phi'*P1*Phi))*beta+Phi'*P1*Phi));
            %P1 = max(eps, P1);
            epsilon(n_step)=(z(n_step)-[thetas{1}(n_step),thetas{2}(n_step),thetas{3}(n_step),thetas{4}(n_step),thetas{5}(n_step),thetas{6}(n_step),thetas{7}(n_step),thetas{8}(n_step),thetas{9}(n_step),thetas{10}(n_step),thetas{11}(n_step)]*Phi)/(1+alpha*(Phi'*Phi)+(Phi'*P1*Phi));    

            for i = 1:11
                %reg_term = lambda * thetas{i}(n_step);
                thetas{i}(n_step+1)=thetas{i}(n_step)+P1(i,:)*Phi*epsilon(n_step+1)-sigma*P1(i,:)*[thetas{1}(n_step),thetas{2}(n_step),thetas{3}(n_step),thetas{4}(n_step),thetas{5}(n_step),thetas{6}(n_step),thetas{7}(n_step),thetas{8}(n_step),thetas{9}(n_step),thetas{10}(n_step),thetas{11}(n_step)]';
                thetas{i}(n_step+1) = max(min(thetas{i}(n_step + 1), UB1(i)), LB1(i));
            end
            X = [thetas{1}(n_step+1), thetas{2}(n_step+1), thetas{3}(n_step+1),thetas{4}(n_step+1), thetas{5}(n_step+1), thetas{6}(n_step+1),thetas{7}(n_step+1), thetas{8}(n_step+1), thetas{9}(n_step+1),thetas{10}(n_step+1),thetas{11}(n_step+1)];
            Ad=[0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0,0;0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1;X(3),X(2),X(6),X(5),X(8),X(7),X(9),X(11),X(10)];
            Bd=[0,0;1,0;0,0;0,1;0,0;0,0;0,0;0,0;X(1),X(4)];
            P=Ad'*P*Ad-(Ad'*P*Bd)*((Bd'*P*Bd+R)\(Bd'*P*Ad))+Q;
            K=(Bd'*P*Bd+R)\(Bd'*P*Ad);
            %K=lqr(Ad,Bd,Q,R);
            %xd(:,n_step)=[delta_U1(max(n_step-2,1)),delta_U1(max(n_step-1,1)),delta_U2(max(n_step-2,1)),delta_U2(max(n_step-1,1)),delta_T0(max(n_step-1,1)),delta_T0(n_step),delta_G(n_step),e(max(n_step-1,1)),e(n_step)]';
            %delta_U1(n_step)=-K(1,:)*xd(:,n_step);
            %delta_U2(n_step)=-K(2,:)*xd(:,n_step);
            %U1(n_step)=U1(max(n_step-1,1))+delta_U1(n_step);
            %U2(n_step)=U2(max(n_step-1,1))+delta_U2(n_step);
            
            xd(:,n_step)=[delta_U1(max(n_step-2,1)),delta_U1(max(n_step-1,1)),delta_U2(max(n_step-2,1)),delta_U2(max(n_step-1,1)),delta_T_amb(max(n_step-1,1)),delta_T_amb(n_step),delta_G_solar(n_step),e(max(n_step-1,1)),e(n_step)]';
            delta_U1(n_step)=-K(1,:)*xd(:,n_step);
            delta_U2(n_step)=-K(2,:)*xd(:,n_step);
            U1(n_step+1)=U1(max(n_step,1))+delta_U1(n_step);
            U2(n_step+1)=U2(max(n_step,1))+delta_U2(n_step);
            %U1(n_step+1)=U1(n_step)+X(1)*e(n_step)+X(2)*(e(n_step)-e(max(n_step-1,1)));%zone
            %U2(n_step+1)=U2(n_step)+X(3)*e(n_step)+X(4)*(e(n_step)-e(max(n_step-1,1)));%tank
        else 
        end
    end 
    total=sqrt(sum(e.^2+0.0002*U1.^2+0.000005*U2.^2)); 
    end
end
