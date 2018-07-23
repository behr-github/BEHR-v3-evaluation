% Make tables of Leitao slopes

load('~/Documents/MATLAB/NO2 Profiles/Workspaces/LeitaoFigures.mat'); 

Coinc = struct('Description',{{}},'All_data_slope',[],'All_data_R2',[],'All_data_std_dev',[],'Coarse_mode_slope',[],'Coarse_mode_R2',[],'Coarse_mode_std_dev',[],'Fine_mode_slope',[],'Fine_mode_R2',[],'Fine_mode_std_dev',[]);
C=1;
AerAbove = Coinc;
A=1;
NO2Above = Coinc;
N=1;

for a=1:6
    for b=1:numel(LeitaoFigures(a).Case)
        if strcmp(LeitaoFigures(a).Case(b).classification,'coincident')
            Coinc.Description{C} = LeitaoFigures(a).Case(b).description;
            Coinc.All_data_slope(C) = LeitaoFigures(a).Case(b).slope;
            Coinc.All_data_R2(C) = LeitaoFigures(a).Case(b).r2;
            Coinc.All_data_std_dev(C) = LeitaoFigures(a).Case(b).sigma_m;
            Coinc.Coarse_mode_slope(C) = LeitaoFigures(a).Case(b).coarse_slope;
            Coinc.Coarse_mode_R2(C) = LeitaoFigures(a).Case(b).coarse_R2;
            Coinc.Coarse_mode_std_dev(C) = LeitaoFigures(a).Case(b).coarse_sigma_m;
            Coinc.Fine_mode_slope(C) = LeitaoFigures(a).Case(b).fine_slope;
            Coinc.Fine_mode_R2(C) = LeitaoFigures(a).Case(b).fine_R2;
            Coinc.Fine_mode_std_dev(C) = LeitaoFigures(a).Case(b).fine_sigma_m;
            C=C+1;
        elseif strcmp(LeitaoFigures(a).Case(b).classification,'aerosol_above');
            AerAbove.Description{A} = LeitaoFigures(a).Case(b).description;
            AerAbove.All_data_slope(A) = LeitaoFigures(a).Case(b).slope;
            AerAbove.All_data_R2(A) = LeitaoFigures(a).Case(b).r2;
            AerAbove.All_data_std_dev(A) = LeitaoFigures(a).Case(b).sigma_m;
            AerAbove.Coarse_mode_slope(A) = LeitaoFigures(a).Case(b).coarse_slope;
            AerAbove.Coarse_mode_R2(A) = LeitaoFigures(a).Case(b).coarse_R2;
            AerAbove.Coarse_mode_std_dev(A) = LeitaoFigures(a).Case(b).coarse_sigma_m;
            AerAbove.Fine_mode_slope(A) = LeitaoFigures(a).Case(b).fine_slope;
            AerAbove.Fine_mode_R2(A) = LeitaoFigures(a).Case(b).fine_R2;
            AerAbove.Fine_mode_std_dev(A) = LeitaoFigures(a).Case(b).fine_sigma_m;
            A=A+1;
        elseif strcmp(LeitaoFigures(a).Case(b).classification,'no2_above');
            NO2Above.Description{N} = LeitaoFigures(a).Case(b).description;
            NO2Above.All_data_slope(N) = LeitaoFigures(a).Case(b).slope;
            NO2Above.All_data_R2(N) = LeitaoFigures(a).Case(b).r2;
            NO2Above.All_data_std_dev(N) = LeitaoFigures(a).Case(b).sigma_m;
            NO2Above.Coarse_mode_slope(N) = LeitaoFigures(a).Case(b).coarse_slope;
            NO2Above.Coarse_mode_R2(N) = LeitaoFigures(a).Case(b).coarse_R2;
            NO2Above.Coarse_mode_std_dev(N) = LeitaoFigures(a).Case(b).coarse_sigma_m;
            NO2Above.Fine_mode_slope(N) = LeitaoFigures(a).Case(b).fine_slope;
            NO2Above.Fine_mode_R2(N) = LeitaoFigures(a).Case(b).fine_R2;
            NO2Above.Fine_mode_std_dev(N) = LeitaoFigures(a).Case(b).fine_sigma_m;
            N=N+1;
        end
    end
end

Coinc.Description = Coinc.Description';
Coinc.All_data_slope = Coinc.All_data_slope';
Coinc.All_data_R2 = Coinc.All_data_R2';
Coinc.All_data_std_dev = Coinc.All_data_std_dev';
Coinc.Coarse_mode_slope = Coinc.Coarse_mode_slope';
Coinc.Coarse_mode_R2 = Coinc.Coarse_mode_R2';
Coinc.Coarse_mode_std_dev = Coinc.Coarse_mode_std_dev';
Coinc.Fine_mode_slope = Coinc.Fine_mode_slope';
Coinc.Fine_mode_R2 = Coinc.Fine_mode_R2';
Coinc.Fine_mode_std_dev = Coinc.Fine_mode_std_dev';

AerAbove.Description = AerAbove.Description';
AerAbove.All_data_slope = AerAbove.All_data_slope';
AerAbove.All_data_R2 = AerAbove.All_data_R2';
AerAbove.All_data_std_dev = AerAbove.All_data_std_dev';
AerAbove.Coarse_mode_slope = AerAbove.Coarse_mode_slope';
AerAbove.Coarse_mode_R2 = AerAbove.Coarse_mode_R2';
AerAbove.Coarse_mode_std_dev = AerAbove.Coarse_mode_std_dev';
AerAbove.Fine_mode_slope = AerAbove.Fine_mode_slope';
AerAbove.Fine_mode_R2 = AerAbove.Fine_mode_R2';
AerAbove.Fine_mode_std_dev = AerAbove.Fine_mode_std_dev';

NO2Above.Description = NO2Above.Description';
NO2Above.All_data_slope = NO2Above.All_data_slope';
NO2Above.All_data_R2 = NO2Above.All_data_R2';
NO2Above.All_data_std_dev = NO2Above.All_data_std_dev';
NO2Above.Coarse_mode_slope = NO2Above.Coarse_mode_slope';
NO2Above.Coarse_mode_R2 = NO2Above.Coarse_mode_R2';
NO2Above.Coarse_mode_std_dev = NO2Above.Coarse_mode_std_dev';
NO2Above.Fine_mode_slope = NO2Above.Fine_mode_slope';
NO2Above.Fine_mode_R2 = NO2Above.Fine_mode_R2';
NO2Above.Fine_mode_std_dev = NO2Above.Fine_mode_std_dev';


LeitaoTables.Coincident = struct2table(Coinc);
LeitaoTables.AerosolAbove = struct2table(AerAbove);
LeitaoTables.NO2Above = struct2table(NO2Above);