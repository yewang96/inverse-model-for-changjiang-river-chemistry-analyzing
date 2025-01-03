clear
clc
tic
data = importdata('major_ions.xlsx'); %load the data
major_ions = data.data;
K = major_ions(:,1);
Ca = major_ions(:,2);
Na = major_ions(:,3);
Mg = major_ions(:,4);
SO4 = major_ions(:,5);
%NO3 = major_ions(:,6);
Re = major_ions(:,7);
sigma =major_ions(:,8); 
runoff = major_ions(:,9);
n = 1000000;

Nsam = length(K);
member = importdata('end_member.xlsx');
end_member = member.data;
Nm = length(end_member);


MemCa_carb = NaN(n,Nsam);
MemMg_carb = NaN(n,Nsam);
MemNa_carb = NaN(n,Nsam);

MemCa_eva = NaN(n,Nsam);
MemMg_eva = NaN(n,Nsam);
MemNa_eva = NaN(n,Nsam);
MemSO4_eva = NaN(n,Nsam);

MemK_ind = NaN(n,Nsam);
MemCa_ind = NaN(n,Nsam);
MemMg_ind = NaN(n,Nsam);
MemNa_ind = NaN(n,Nsam);
MemSO4_ind = NaN(n,Nsam);
%MemNO3_ind = NaN(n,Nsam);
MemRe_ind = NaN(n,Nsam);

MemK_farm = NaN(n,Nsam);
MemCa_farm = NaN(n,Nsam);
MemMg_farm = NaN(n,Nsam);
MemNa_farm = NaN(n,Nsam);
MemSO4_farm = NaN(n,Nsam);
%MemNO3_farm = NaN(n,Nsam);

MemK_shale = NaN(n,Nsam);
MemCa_shale = NaN(n,Nsam);
MemMg_shale = NaN(n,Nsam);
MemNa_shale = NaN(n,Nsam);
MemSO4_shale = NaN(n,Nsam);
%MemNO3_shale = NaN(n,Nsam);
MemRe_shale = NaN(n,Nsam);

MemK_coal = NaN(n,Nsam);
MemCa_coal = NaN(n,Nsam);
MemMg_coal = NaN(n,Nsam);
MemNa_coal = NaN(n,Nsam);
MemSO4_coal = NaN(n,Nsam);
%MemNO3_coal = NaN(n,Nsam);
MemRe_coal = NaN(n,Nsam);

F_carb = NaN(n,length(K));
F_eva = NaN(n,length(K));
F_ind = NaN(n,length(K));
F_farm = NaN(n,length(K));
F_coal = NaN(n,length(K));
F_shale = NaN(n,length(K));


for i = 1:Nm
    pd = makedist('Uniform','lower',end_member(i,1),'upper',end_member(i,2));
    eval(['pd_',member.textdata{i+1,1}, '= pd;']);
end

for i = 1:n
    for j = 1:Nm
        eval(['pdi = pd_',member.textdata{j+1,1},';']); 
        eval([member.textdata{j+1,1}, '= random(pdi);']); 
    end
    Mem = [ 0,       0,        K_ind,      K_coal,    K_shale; 
            Ca_eva,  Ca_carb,  Ca_ind,     Ca_coal,   Ca_shale;
            Na_eva,  Na_carb,  Na_ind,     Na_coal,   Na_shale;
            Mg_eva,  Mg_carb,  Mg_ind,     Mg_coal,   Mg_shale;
            SO4_eva, 0,        SO4_ind,    SO4_coal,  SO4_shale;
            0,       0,        Re_ind,     Re_coal,   Re_shale ];
         
    Mar = [K';Ca';Na';Mg';SO4';Re'];
    F = Mem\Mar;

    %verification
    neg = find(F < 0);
    F(neg) = 1000;
    F_sum = sum(F,1); %sum absolute values of each row
	F_indi = find(F_sum < 2); %find indices where sum(x)=sum(abs(x))
        

    F_eva(i,F_indi) = F(1,F_indi); %mixing fractions     
    F_carb(i,F_indi) = F(2,F_indi);
    F_ind(i,F_indi) = F(3,F_indi);
    F_coal(i,F_indi) = F(4,F_indi);
    F_shale(i,F_indi) = F(5,F_indi);
    
    Class = ~isnan(F_eva(i,:));
    for pp = 1:length(Class)
        if Class(pp) == 1
           MemCa_eva(i,pp) = Mem(2,1);
           MemNa_eva(i,pp) = Mem(3,1);
           MemMg_eva(i,pp) = Mem(4,1);
           MemSO4_eva(i,pp) = Mem(5,1);

           MemCa_carb(i,pp) = Mem(2,2);
           MemNa_carb(i,pp) = Mem(3,2);
           MemMg_carb(i,pp) = Mem(4,2);

           MemK_ind(i,pp) = Mem(1,3);
           MemCa_ind(i,pp) = Mem(2,3);
           MemNa_ind(i,pp) = Mem(3,3);
           MemMg_ind(i,pp) = Mem(4,3);
           MemSO4_ind(i,pp) = Mem(5,3);
           %MemNO3_ind(i,pp) = Mem(6,3);
           MemRe_ind(i,pp) = Mem(6,3);
           
           MemK_coal(i,pp) = Mem(1,4);
           MemCa_coal(i,pp) = Mem(2,4);
           MemNa_coal(i,pp) = Mem(3,4);
           MemMg_coal(i,pp) = Mem(4,4);
           MemSO4_coal(i,pp) = Mem(5,4);
           %MemNO3_coal(i,pp) = Mem(6,5);
           MemRe_coal(i,pp) = Mem(6,4);

           MemK_shale(i,pp) = Mem(1,5);
           MemCa_shale(i,pp) = Mem(2,5);
           MemNa_shale(i,pp) = Mem(3,5);
           MemMg_shale(i,pp) = Mem(4,5);
           MemSO4_shale(i,pp) = Mem(5,5);
           %MemNO3_shale(i,pp) = Mem(6,6);
           MemRe_shale(i,pp) = Mem(6,5);


        end
    end
    
end

F_eva_mean(1,:) = mean(F_eva,"omitnan");
F_carb_mean(1,:) = mean(F_carb,"omitnan");
F_ind_mean(1,:) = mean(F_ind,"omitnan");
F_coal_mean(1,:) = mean(F_coal,"omitnan");
F_shale_mean(1,:) = mean(F_shale,"omitnan");



%% CO2 sink through carbonate weathering
cation_carb = F_carb.*MemCa_carb + F_carb.*MemMg_carb ; 

cation_carb_cc = cation_carb'.*sigma.*runoff;    % (mol/km2/yr)
cation_carb_p = prctile(cation_carb_cc,[25,50,75],2);  

carb_weather = (F_carb.*MemCa_carb*20 + F_carb.*MemMg_carb*12)'.*sigma.*runoff/1000000;
carb_weather_p = prctile(carb_weather,[25,50,75],2);  % 25th,50th,75th value of data (t/km2/yr)


%% Re source portioning
Re_ind_cs = F_ind.*MemRe_ind;  % eq/eq
Re_coal_cs = F_coal.*MemRe_coal;  % eq/eq
Re_shale_cs = F_shale.*MemRe_shale;  % eq/eq

Re_ind_conc = Re_ind_cs'.*sigma;  % ( pmol/L)
Re_ind_conp = prctile(Re_ind_conc,[25,50,75],2);
Re_ind_cc = Re_ind_cs'.*sigma.*runoff;  % ( umol/km2/yr)
Re_ind_p = prctile(Re_ind_cc,[25,50,75],2);   % 25th,50th,75th value of data (J_Re_ind)

Re_coal_conc = Re_coal_cs'.*sigma;  % ( pmol/L)
Re_coal_conp = prctile(Re_coal_conc,[25,50,75],2);
Re_coal_cc = Re_coal_cs'.*sigma.*runoff;  % ( umol/km2/yr)
Re_coal_p = prctile(Re_coal_cc,[25,50,75],2);   % 25th,50th,75th value of data (J_Re_coal)

Re_shale_conc = Re_shale_cs'.*sigma;  % ( pmol/L)
Re_shale_conp = prctile(Re_shale_conc,[25,50,75],2);
Re_shale_cc = Re_shale_cs'.*sigma.*runoff;  % ( umol/km2/yr)
Re_shale_p = prctile(Re_shale_cc,[25,50,75],2);   % 25th,50th,75th value of data (J_Re_shale)
Re_shale_mean = mean(Re_shale_conc,"omitnan");

%% CO2 yield
J_shale = NaN(n,Nsam);
pro = 1;
for col = 1:Nsam
    
    F_gra = 0.7 + (1-0.7)*rand(n,1);
    RevsOC = 3.16 + (5.65-3.16)*rand(n,1);


    J_shale_r = Re_shale_cc(col,:).*(185/100000)./RevsOC'.*F_gra';  %tC/km2/yr
    J_shale(:,col) = J_shale_r';  %tC/km2/yr

    if sum(isnan(J_shale_r)) ~= n
       figure (1)
       subplot(5,5,pro),histogram(J_shale_r,100,'FaceColor','r');  %normal distribuion test
       xlabel('OCpetro oxidation rate (tC/km2/yr)');
       ylabel('Frequency');
       hold on;
       pro = pro+1;
  
    end

end 
J_shale_mean =  mean(J_shale,"omitnan");   %mean value
J_shale_std = std(J_shale,"omitnan");  %std value
J_shale_p = prctile(J_shale',[25,50,75],2);   % 25th,50th,75th value of data  (tC/km2/yr)



toc




