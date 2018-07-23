% Run_all_aer_categories
%   Simple script that will run Run_Spiral_Verification for each of the six
%   categories regarding relative distribution of aerosol and NO2
%
%   Josh Laughner <joshlaugh5@gmail.com> 30 Jan 2015

% A reminder to the user, and a chance to cancel 
msg = 'MAKE SURE RUN_SPIRAL_VERIFICATION IS SET TO USE THE CORRECT CAMPAIGN!!!';
msg2 = 'Waiting 4 sec - Press Ctrl+C to cancel';
l = length(msg); l2 = length(msg2); lspc = floor((l-l2)/2);
spc = repmat(' ',1,lspc);
msg2 = sprintf('%s%s%s',spc,msg2,spc);
accent = repmat('*',1,length(msg));
fprintf('%s\n%s\n%s\n%s\n',accent,msg,msg2,accent);
pause(4);

% Make this equal to your structure with the category information.
CatStruct = DAQ_CO_cat;

% The six categories to run
aer_cats = {'CoincidentLow','CoincidentHigh','AerosolAboveLow','AerosolAboveHigh','NO2AboveLow','NO2AboveHigh'};

% Prepare a 1x6 structure to receive the outputs
Comparison = struct('Category',aer_cats,'lon_iall',cell(1,6),'lat_iall',cell(1,6),...
    'airno2_iall',cell(1,6),'omino2_iall',cell(1,6),'behrno2_iall',cell(1,6),...
    'db_iall',struct,'dates_iall',cell(1,6));

% Loop through the six categories, calling Run_Spiral_Verification for each
for c=1:6
    [Comparison(c).lon_iall, Comparison(c).lat_iall, Comparison(c).omino2_iall, Comparison(c).behrno2_iall,...
        Comparison(c).airno2_iall, Comparison(c).db_iall, Comparison(c).dates_iall] =...
        Run_Spiral_Verification('profnums', CatStruct.(aer_cats{c}));
end