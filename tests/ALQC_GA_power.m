tic 
clear     
%weather information (Weather_and_ElectrictyPrice_DATA.xlsx) 
%year 2015 in Finland
T_amb_temp=readmatrix('Helsinki_reference_year_new.xlsx','Range','B1730:B1826');%Hourly outdoor temperature on the seventy-fifth day
G_temp=readmatrix('Helsinki_reference_year_new.xlsx','Range','H1730:H1826');%Hourly solar radiation on the seventy-fifth day
alpha_int=0.5;alpha_sol=0.3;alpha_inf=0.8;
Num_B1=33;Num_B2=33;Num_B3=34;%number of B1,B2,B3
DHW_WD_March=43*1.020;%Average domestic heat and water consumption per person in March（kg
hourlyfactor_fix=[0,0,0,0,0,0,0,0.275,0.065,0.04,0.02,0.04,0.05,0,0.02,0.02,0.02,0,0.05,0.02,0.13,0.26,0,0];%every day
% load('data1_Gather_P_norm_5min.mat','G','To');
load("random_value1.mat","flag_B1","Occupants","factor_LW","flag_B2","Occupants_EL","factor_EL","flag_B3","Occupants_MP","factor_MP","T_g","ini_zone_B1","ini_tank_B1","T_g_EL","ini_zone_B2","ini_tank_B2","T_g_MP","ini_zone_B3","ini_tank_B3");
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%building 1:LW building%%%%%%%%%%%%%%%%%%%%%%%
%flag_B1=0.6*rand(1,Num_B1)+0.7;
aa=8.5*flag_B1;%m Num_B1
bb=10.6*flag_B1;%m
A_roof=aa.*bb;%m^2
A_ground=A_roof;%1*Num_B1


%internal heat gain
%Occupants=randi([2,4],1,Num_B1);%Randomly generate a value of 2, 3 and 4, which represents the number of people living in the building
Q_people_day=zeros(24,Num_B1);%internal heat gain of people in weekday
Q_light_day=zeros(24,Num_B1);%internal heat gain of light in weekday
Q_equipment_day=zeros(24,Num_B1);%internal heat gain of equipment in weekday
for j=1:Num_B1
    for i=1:7
        Q_people_day(i,j)=(126*Occupants(j));
        Q_equipment_day(i,j)=15*2*A_ground(j)*0.1;
    end
    Q_people_day(8,j)=(126*Occupants(j))*0.2;
    Q_equipment_day(8,j)=15*2*A_ground(j)*0.25;
    for i=9:15
         Q_equipment_day(i,j)=15*2*A_ground(j)*0.1;
    end
    for i=16:17
        Q_people_day(i,j)=(126*Occupants(j))*0.25;
        Q_equipment_day(i,j)=15*2*A_ground(j)*0.2;
    end
    for i=16:18
        Q_light_day(i,j)=10*2*A_ground(j)*0.1;
    end
    for i=19:23
        Q_light_day(i,j)=10*2*A_ground(j)*0.5;
    end
    for i=18:23
        Q_people_day(i,j)=(126*Occupants(j));
        Q_equipment_day(i,j)=15*2*A_ground(j)*0.5;
    end
    Q_people_day(24,j)=(126*Occupants(j));
    Q_equipment_day(24,j)=15*2*A_ground(j)*0.1;
    
end
Q_int_day=Q_people_day+Q_light_day*alpha_int+Q_equipment_day*alpha_int;

%%%domestic hot water load
hourlyfactor_LW=ones(Num_B1,1)*hourlyfactor_fix;%Generate a matrix with hourlyfactor_fix vector in each row, and the number of rows is the total number of buildings
%factor_LW=0.004*rand(Num_B1,24)-0.002;
hourlyfactor_LW=hourlyfactor_LW+factor_LW;%give random factor to DHW (Num*24 D)
m1_LW=zeros(Num_B1,24);
for i=1:Num_B1
    for j=1:24
    if hourlyfactor_LW(i,j)<0
       hourlyfactor_LW(i,j)=0;
    else
    end
    m1_LW(i,j)=DHW_WD_March*Occupants(i)*hourlyfactor_LW(i,j);
    end
end

%%  %%%%%%%%%%%%%%%%%%%%%%%building 2:EL building%%%%%%%%%%%%%%%%%%%%%%%
%flag_B2=0.6*rand(1,Num_B2)+0.7;
cc=15*flag_B2;dd=10*flag_B2;
A_roof_EL=cc.*dd;%m^2
A_ground_EL=A_roof_EL;
%internal heat gain
%Occupants_EL=randi([2,4],1,Num_B2);%Randomly generate a value of 2, 3 and 4, which represents the number of people living in the building
Q_people_day_EL=zeros(24,Num_B2);%internal heat gain of people in weekday
Q_light_day_EL=zeros(24,Num_B2);%internal heat gain of light in weekday
Q_equipment_day_EL=zeros(24,Num_B2);%internal heat gain of equipment in weekday

for j=1:Num_B2
    for i=1:7
        Q_people_day_EL(i,j)=(126*Occupants_EL(j));
        Q_equipment_day_EL(i,j)=15*A_ground_EL(j)*0.1;
    end
    Q_people_day_EL(8,j)=(126*Occupants_EL(j))*0.5;
    Q_equipment_day_EL(8,j)=15*A_ground_EL(j)*0.2;
    for i=9:15
         Q_equipment_day_EL(i,j)=15*A_ground_EL(j)*0.1;
    end
    for i=16:17
        Q_people_day_EL(i,j)=(126*Occupants_EL(j))*0.25;
        Q_equipment_day_EL(i,j)=15*A_ground_EL(j)*0.5;
    end
    for i=16:18
        Q_light_day_EL(i,j)=10*A_ground_EL(j)*0.25;
    end
    for i=19:23
        Q_light_day_EL(i,j)=10*A_ground_EL(j)*0.5;
    end
    for i=18:23
        Q_people_day_EL(i,j)=(126*Occupants_EL(j));
        Q_equipment_day_EL(i,j)=15*A_ground_EL(j)*0.5;
    end
    Q_people_day_EL(24,j)=(126*Occupants_EL(j));
    Q_equipment_day_EL(24,j)=15*A_ground_EL(j)*0.1;
end
Q_int_day_EL=Q_people_day_EL+Q_light_day_EL*alpha_int+Q_equipment_day_EL*alpha_int;

%%%domestic hot water load
hourlyfactor_EL=ones(Num_B2,1)*hourlyfactor_fix;%Generate a matrix with hourlyfactor_fix vector in each row, and the number of rows is the total number of buildings
%factor_EL=0.004*rand(Num_B2,24)-0.002;
hourlyfactor_EL=hourlyfactor_EL+factor_EL;%give random factor to DHW (Num*24 D)
m1_EL=zeros(Num_B2,24);
for i=1:Num_B2
    for j=1:24
    if hourlyfactor_EL(i,j)<0
       hourlyfactor_EL(i,j)=0;
    else
    end
    m1_EL(i,j)=DHW_WD_March*Occupants_EL(i)*hourlyfactor_EL(i,j);
    end
end

%%  %%%%%%%%%%%%%%%%%building 3:MP building%%%%%%%%%%%%%%%%%%%%%%%
%flag_B3=0.6*rand(1,Num_B3)+0.7;
ee=8.5*flag_B3;ff=10.6*flag_B3;
A_roof_MP=ee.*ff;%m^2
A_ground_MP=A_roof_MP;
%internal heat gain
%Occupants_MP=randi([2,4],1,Num_B3);%Randomly generate a value of 2, 3 and 4, which represents the number of people living in the building
Q_people_day_MP=zeros(24,Num_B3);%internal heat gain of people in weekday
Q_light_day_MP=zeros(24,Num_B3);%internal heat gain of light in weekday
Q_equipment_day_MP=zeros(24,Num_B3);%internal heat gain of equipment in weekday

for j=1:Num_B3
    for i=1:7
        Q_people_day_MP(i,j)=(126*Occupants_MP(j));
        Q_equipment_day_MP(i,j)=15*2*A_ground_MP(j)*0.1;
    end
    Q_people_day_MP(8,j)=(126*Occupants_MP(j))*0.5;
    Q_equipment_day_MP(8,j)=15*2*A_ground_MP(j)*0.2;
    for i=9:15
         Q_equipment_day_MP(i,j)=15*2*A_ground_MP(j)*0.1;
    end
    for i=16:17
        Q_people_day_MP(i,j)=(126*Occupants_MP(j))*0.25;
        Q_equipment_day_MP(i,j)=15*2*A_ground_MP(j)*0.2;
    end
    for i=16:18
        Q_light_day_MP(i,j)=10*2*A_ground_MP(j)*0.1;
    end
    for i=19:23
        Q_light_day_MP(i,j)=10*2*A_ground_MP(j)*0.5;
    end
    for i=18:23
        Q_people_day_MP(i,j)=(126*Occupants_MP(j));
        Q_equipment_day_MP(i,j)=15*2*A_ground_MP(j)*0.5;
    end
    Q_people_day_MP(24,j)=(126*Occupants_MP(j));
    Q_equipment_day_MP(24,j)=15*2*A_ground_MP(j)*0.1;

end
Q_int_day_MP=Q_people_day_MP+Q_light_day_MP*alpha_int+Q_equipment_day_MP*alpha_int;

%%%domestic hot water load
hourlyfactor_MP=ones(Num_B3,1)*hourlyfactor_fix;%Generate a matrix with hourlyfactor_fix vector in each row, and the number of rows is the total number of buildings
%factor_MP=0.004*rand(Num_B3,24)-0.002;
hourlyfactor_MP=hourlyfactor_MP+factor_MP;%give random factor to DHW (Num*24 D)
m1_MP=zeros(Num_B3,24);
for i=1:Num_B3
    for j=1:24
    if hourlyfactor_MP(i,j)<0
       hourlyfactor_MP(i,j)=0;
    else
    end
    m1_MP(i,j)=DHW_WD_March*Occupants_MP(i)*hourlyfactor_MP(i,j);
    end
end

%%  %%%%%%%%%%%%%%%%%  GA  %%%%%%%%%%%%%%%%%%%%%%%

options = optimoptions('ga','PopulationSize',250);
LB= [-1 -1 -6 2 0 3 1 1];
UB= [0 0 -2 5 3 6 4 4];
% UB=[10*10^(-2),10*10^(-3),5*10^(-2),10*10^(-2),10*10^(-3),5*10^(-2)];
[theta,total] = ga(@(A)ALQC_GA_function_power(A,G_temp,T_amb_temp,Q_int_day,Q_int_day_EL,Q_int_day_MP,m1_LW,m1_EL,m1_MP,flag_B1,flag_B2,flag_B3,Num_B1,Num_B2,Num_B3,T_g,T_g_EL,T_g_MP,ini_zone_B1,ini_zone_B2,ini_zone_B3,ini_tank_B1,ini_tank_B2,ini_tank_B3),8,[],[],[],[],LB,UB,[],options);
para=theta;
P_total=total;
save('ALQC_GA_power.mat','para','P_total');

% load('PID_GA_para.mat','para','P_total');
% figure(2)
% plot(y','linewidth',1);
% hold on
% plot(yd','linewidth',1);
% xlabel('Time (h)');
% ylabel('Normalized average aggregated power');
% legend('z_f','y_d');
% xlabel('Time (Hour)');
% axis([0 288*4 0 0.6]);
% set(gca,'xtick',0:12*4:288*4 )
% set(gca,'Xticklabel',{'0','4','8','12','16','20','24','28','32','36','40','44','48','52','56','60','64','68','72','76','80','84','88','92','96'});
% 
% 
% figure(3)
% plot(U1,'linewidth',1);
% hold on
% plot(U2,'linewidth',1);
% xlabel('Time (h)');
% ylabel('Input (℃)');
% axis([0 288*4 -5 5]);
% legend('Input of the zone U_1','Input of the tank U_2');
% set(gca,'xtick',0:12*4:288*4 )
% set(gca,'Xticklabel',{'0','4','8','12','16','20','24','28','32','36','40','44','48','52','56','60','64','68','72','76','80','84','88','92','96'});

toc