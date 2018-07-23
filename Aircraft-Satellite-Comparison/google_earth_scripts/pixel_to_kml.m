%pixcorn_kml: Write pixel boundaries to a KML file
latcorn = db.latcorn;
loncorn = db.loncorn;

folder = '/Users/Josh/Desktop/kml_test';
filenames = {'Pixels 7-1-2011.kml'};
colors = {'black','white','black','white','black','white','black','white'};

if iscell(latcorn)
    n=numel(latcorn);
    for a=1:n
        kml = geoshape();
        m = size(latcorn{a},2);
        for b=1:m
            kml(b).Latitude = [latcorn{a}(:,b);latcorn{a}(1,b)];
            kml(b).Longitude = [loncorn{a}(:,b);loncorn{a}(1,b)];
            kml(b).Color = colors{a};
        end
        if ~isempty(kml)
            fname = fullfile(folder,filenames{a});
            kmlwrite(fname,kml,'Color',kml.Color,'Width',4);
        end
    end
else
    kml = geoshape();
    m = size(latcorn,2);
    for b=1:m
        kml(b).Latitude = [latcorn(:,b);latcorn(1,b)];
        kml(b).Longitude = [loncorn(:,b);loncorn(1,b)];
        kml(b).Color = colors{1};
    end
    if ~isempty(kml)
        fname = fullfile(folder,filenames{1});
        kmlwrite(fname,kml,'Color',kml.Color,'Width',4);
    end
end