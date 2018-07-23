% Calculate the coarse and fine only slopes for LeitaoFigures

load('~/Documents/MATLAB/NO2 Profiles/Workspaces/LeitaoFigures.mat'); % loads the LeitaoFigures variable

for a=1:6
    for b=1:numel(LeitaoFigures(a).Case)
        amfs = LeitaoFigures(a).Case(b).amf_over_amf0_percent(LeitaoFigures(a).Case(b).coarse_mode);
        aods = LeitaoFigures(a).Case(b).aod(LeitaoFigures(a).Case(b).coarse_mode);
        if ~isempty(amfs)
            [~,~,~,D] = calc_fit_line(aods,amfs,'regression','rma');
            LeitaoFigures(a).Case(b).coarse_slope = D.P(1);
            LeitaoFigures(a).Case(b).coarse_R2 = D.R2;
            LeitaoFigures(a).Case(b).coarse_sigma_m = D.StdDevM;
        else
            LeitaoFigures(a).Case(b).coarse_slope = nan;
            LeitaoFigures(a).Case(b).coarse_R2 = nan;
            LeitaoFigures(a).Case(b).coarse_sigma_m = nan;
        end
        
        amfs = LeitaoFigures(a).Case(b).amf_over_amf0_percent(LeitaoFigures(a).Case(b).fine_mode);
        aods = LeitaoFigures(a).Case(b).aod(LeitaoFigures(a).Case(b).fine_mode);
        if ~isempty(amfs)
            [~,~,~,D] = calc_fit_line(aods,amfs,'regression','rma');
            LeitaoFigures(a).Case(b).fine_slope = D.P(1);
            LeitaoFigures(a).Case(b).fine_R2 = D.R2;
            LeitaoFigures(a).Case(b).fine_sigma_m = D.StdDevM;
        else
            LeitaoFigures(a).Case(b).fine_slope = nan;
            LeitaoFigures(a).Case(b).fine_R2 = nan;
            LeitaoFigures(a).Case(b).fine_sigma_m = nan;
        end
    end
end