clc 
clear
len=2; %distance in km.
bred=2; %distance in km
Area=len*bred;
digits(8);
Data = readtable('Water_Data3.csv','ReadRowNames',true,'Format','%u%f%f');
S = size(Data,1); %Number of sensors
Sensor = zeros(S,2);
disp('Number of Sensors deployed:')
disp(S)
Data_BS = readtable('BS_Locations3.txt','ReadRowNames',true,'Format','%u%f%f');
B = size(Data_BS,1); %Number of BSs
BaseStation = zeros(B,2);
rad = 0.65; %This coverage radius of BS in km, corresponds to 10^-3 (0.001) outage prob. threshold.
disp('Total Number of BSs used:')
disp(B)
disp('Coverage radius of each BS:')
disp(rad)

for s=1:S
    Sensor(s,1) = vpa(Data.x(s)); %vpa function, along with digits is used to extract precision 6 data from table.
    Sensor(s,2) = vpa(Data.y(s));
end

Coverage_BS = zeros(1,B);
for b=1:B
    BaseStation(b,1) = vpa(Data_BS.x(b));
    BaseStation(b,2) =  vpa(Data_BS.y(b));   
    Coverage_BS(b) = rad; %pi*rad^2;
end

X_b_s = zeros(S,B);
for s=1:S
    for b=1:B       
        dist = Long_Lat_Dist(Sensor(s,:),BaseStation(b,:));
        if (dist <= Coverage_BS(b))
            X_b_s(s,b) = 1;
        end
    end
end

%Check for constraint satisfaction first
Not_Covered =zeros(0,2);
m=0;
for s=1:S
    if(sum(X_b_s(s,:)) >= 1)
        continue
    else
        m=m+1;
        Not_Covered(m,:) = Sensor(s,:);      
    end
end

if (~(isempty(Not_Covered)))
    disp('Not Sufficient Base stations deployed')
    return
else
    
    %Greedy Algo
    BS_Deploy = zeros(0,2);
    Z_b_s = X_b_s;
    while (sum(Z_b_s(:))) %Iterate untill all sensors are covered
        [Max_Val,Max_Idx] = max(sum(Z_b_s,1)); %sum(Z_b_s,1) gives number of sensors covered by each BS
        BS_Deploy = [BS_Deploy;BaseStation(Max_Idx,:)]; %Appending the BS which covers the most sensors at this iteration.
        for s=1:S
            if (Z_b_s(s,Max_Idx) == 1)
                Z_b_s(s,:) = zeros(1,B); %Removing the sensors already covered the selected BS
            end
        end
    end
   disp('Greedy Algo. BS count:') 
   disp(size(BS_Deploy,1)) 
   disp('Greedy Algo. BS locations:')
   disp(BS_Deploy)
   
    %Exhaustive Search Loop
    BS_Vector = 1:B;
    for k=1:B    
        BS_Combo = nchoosek(BS_Vector,k); %Matrix of all possible combinations https://www.mathworks.com/help/matlab/ref/nchoosek.html
        for i=1:nchoosek(B,k) %Number of all possible combinations
            Temp_Vector = zeros(S,1);
            for j=1:k
                Temp_Vector = Temp_Vector | X_b_s(:,BS_Combo(i,j)); %OR operation btw Base stations to collect sensors covered by them
            end
            Sensor_Count = sum(Temp_Vector); %Number of sensors covered by the current BS combination
            if (Sensor_Count == S)
                disp('Optimal BS count:')
                disp(k)
                disp('Optimal BS locations:')                
                disp(BaseStation(BS_Combo(i,:),:))
                figure()
                grid on
                hold on 
                scatter(BaseStation(:,1),BaseStation(:,2),'x','LineWidth',1)
                scatter(Sensor(:,1),Sensor(:,2),'LineWidth',0.5)
                scatter(BaseStation(BS_Combo(i,:),1),BaseStation(BS_Combo(i,:),2),'o','LineWidth',12,'MarkerEdgeColor','r')
                scatter(BS_Deploy(:,1),BS_Deploy(:,2),'o','LineWidth',6,'MarkerEdgeColor','c')
                hold off                
                set(gca,'FontSize',14)                 
                xlabel('Longitude in degrees');
                ylabel('Latitude in degrees');                
                legend('Base Stations', 'Sensors', 'Optimal', 'Greedy')
                legend('NumColumns',4)
                return
            end
        end    
    end
end % End for if (~(isempty(Not_Covered)))