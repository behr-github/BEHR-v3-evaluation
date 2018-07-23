% Script for the figures of Leitao (2010) profiles from the data in the
% .mat file Leitao 16 Scenarios (saved in NO2 Profiles/Workspaces)
function PlotLeitaoScenarios

%% Make a 8 panel figure for each of the 8 box plots
scenarios = {'A','B','C','D','E','F','G','H'};
descriptor = repmat({''},1,8);
figure;
MakeSubplots(2,4,scenarios,descriptor);

%% Make a three panel figure for the rural and 2 urban plots
scenarios = {'I','J','K'};
descriptor = {' - Rural profiles',' - Urban profiles',' - Urban profiles'};
figure;
MakeSubplots(1,3,scenarios,descriptor);

%% Make a three and two panel figure for the desert dust and biomass burning profiles
scenarios = {'L','M','N'};
descriptor = {' - Desert dust',' - Desert dust',' - Desert dust'};
figure; 
MakeSubplots(1,3,scenarios,descriptor);

scenarios = {'O','P'};
descriptor = {' - Biomass burning',' - Biomass burning'};
figure;
MakeSubplots(1,2,scenarios,descriptor);

end

function MakeSubplots(rows,cols,scenarios,descriptor)
    for a=1:numel(scenarios)
        load('~/Documents/MATLAB/NO2 Profiles/Workspaces/Leitao 16 scenarios.mat'); % introduces the variable LeitaoScenarios
        no2 = LeitaoScenarios.(scenarios{a}).NO2; %
        aer = LeitaoScenarios.(scenarios{a}).AerosolExt;
        alts = LeitaoScenarios.(scenarios{a}).alts; % in km
        % We'll normalized the NO2 and aerosol fields by the total area under
        % the curve before plotting (so we're plotting the shape factor)
        int_no2 = trapz(alts,no2);
        int_aer = trapz(alts,aer);
        no2 = no2 / int_no2;
        aer = aer / int_aer;
        
        subplot(rows,cols,a);
        set(gca,'fontsize',14);
        lno2 = line(no2,alts,'color','b','linestyle','--','linewidth',2);
        laer = line(aer,alts,'color','r','linestyle','-.','linewidth',1);
        
        legend([lno2;laer],{'NO_2','Aerosol'});
        xlabel('Shape factor (conc./column)')
        ylabel('Altitude (km)');
        title(sprintf('Scenario %s%s',scenarios{a},descriptor{a}));
    end
end