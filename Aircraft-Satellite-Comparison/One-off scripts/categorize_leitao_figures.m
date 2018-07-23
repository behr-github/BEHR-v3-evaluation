function [ LeitaoFigures ] = categorize_leitao_figures( LeitaoFigures )
%categorize_leitao_figures Categorize the aerosol/no2 profiles from Leitao et. al. 2010

for a=1:numel(LeitaoFigures)
    for b=1:numel(LeitaoFigures(a).Case)
        category = multiple_categorize_profile(LeitaoFigures(a).Case(b).no2_prof, LeitaoFigures(a).Case(b).prof_alts, LeitaoFigures(a).Case(b).aer_prof, LeitaoFigures(a).Case(b).prof_alts);
        LeitaoFigures(a).Case(b).classification = regexprep(category{1},' ','_');
    end
end

end

