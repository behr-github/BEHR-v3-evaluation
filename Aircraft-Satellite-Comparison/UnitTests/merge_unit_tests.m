function [  ] = merge_unit_tests(  )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

test_combine_merges;

    function status = test_combine_merges()
        [TotMerge, M1M2, M1M3, M1M4, M1, M2, M3, M4] = create_pseudo_merges();
        
        M1M2_test = combine_merges(M1,M2);
        M1M3_test = combine_merges(M1,M3);
        M1M4_test = combine_merges(M1,M4);
        Tot_test = combine_merges(M1,M2,M3,M4);
        
        test_vec = [isequal(M1M2, M1M2_test), isequal(M1M3, M1M3_test), isequal(M1M4, M1M4_test), isequal(TotMerge, Tot_test)];
        test_descr = {'M1 + M2 (different fields, same UTC',...
                      'M1 + M3 (same fields, different UTC',...
                      'M1 + M4 (different fields, different UTC',...
                      'M1 + M2 + M3 + M4 (adding multiple merges at once)'};
        status = all(test_vec);
        fprintf('TEST_COMBINE_MERGES, checking merge returned successfully.\n');
        if ~status
            fprintf('\tFAILED. Causes:\n  %s\n', strjoin(test_descr(~test_vec), '\n\t  '));
        else
            fprintf('\tPASS\n');
        end
        
        fprintf('\nTEST_COMBINE_MERGES, testing error checking...\n');
        M2bad = M2;
        M2bad.Data.UTC.Unit = 'not unitless';
        status = status && error_check_combine_merges(M1, M2bad, 'Merges with different UTC units');
        
        M2bad = M2;
        M2bad.Data.UTC.Fill = -M2bad.Data.UTC.Fill;
        status = status && error_check_combine_merges(M1, M2bad, 'Merges with different UTC fill values');
        
        M2bad = M2;
        M2bad.metadata.date = '2001-05-09'; % random different date
        status = status && error_check_combine_merges(M1, M2bad, 'Merges with different dates');
        
        M2bad = M2;
        M2bad.metadata.upper_lod_flag = -M2bad.metadata.upper_lod_flag;
        status = status && error_check_combine_merges(M1, M2bad, 'Merges with different upper_lod_flags');
        
        M3bad = M3;
        M3bad.Data.alpha.Fill = -M3bad.Data.alpha.Fill;
        status = status && error_check_combine_merges(M1, M3bad, 'Merges with different field fill values');
        
        M3bad = M3;
        M3bad.Data.alpha.Unit = 'not a unit';
        status = status && error_check_combine_merges(M1, M3bad, 'Merges with different field units');
        
        M3bad = M3;
        M3bad.Data.UTC.Values = M1.Data.UTC.Values;
        status = status && error_check_combine_merges(M1, M3bad, 'Merges with values for the same UTC in the same fields');
        
        if status
            stat_str = 'PASS';
        else
            stat_str = 'FAIL';
        end
        fprintf('\nTEST_COMBINE_MERGES: overall %s\n', stat_str);
    end

end

function found_error = error_check_combine_merges(M1, M2, test_descr)
fprintf('\t%s: ', test_descr);
found_error = false;
try
    combine_merges(M1,M2);
catch err
    if strcmp(err.identifier, 'combine_merges:combine_merge_fail')
        found_error = true;
    else
        fprintf('FAIL (did not result in "combine_merge_fail" error)\n');
    end
end
if ~found_error
    fprintf('FAIL (did not result in any error)\n');
else
    fprintf('PASS (gave message "%s")\n', err.message);
end
end

function [TotalMerge, M1M2, M1M3, M1M4, M1, M2, M3, M4] = create_pseudo_merges()
% Create fake merges by making one big merge and breaking it down into four
% merges. This will let us test combining merges with the same UTCs but
% different fields, similar fields but different UTCs, and different fields
% and different UTCs.

mdata = struct('file', '',...
               'description', '',...
               'date', '2000-01-01',...
               'info', '',...
               'upper_lod_flag', -777777,...
               'upper_lod_value', 'N/A',...
               'lower_lod_flag', -888888,...
               'lower_lod_value', 'N/A');
utcs = 1:60;
n = numel(utcs);

field_fill = -9999;
dstr = struct('Unit', 'unitless', 'Fill', field_fill, 'Values', field_fill*ones(1,n));
           
alpha_vals = rand(1,n)*1000;
beta_vals = randi(20,1,n);
gamma_vals = rand(1,n)*256;
delta_vals = randi(256,1,n);



M1.metadata = mdata;
M1.metadata.file = 'tfile1.ict';
M1.metadata.description = 'Testing merge structure 1';
M1.metadata.info = 'UTCs 1-30, fields alpha & beta';

m1_inds = 1:30;
M1.Data.UTC = dstr;
M1.Data.UTC.Values = utcs(m1_inds);
M1.Data.alpha = dstr;
M1.Data.alpha.Values = alpha_vals(m1_inds);
M1.Data.beta = dstr;
M1.Data.beta.Values = beta_vals(m1_inds);

M2.metadata = mdata;
M2.metadata.file = 'tfile2.ict';
M2.metadata.description = 'Testing merge structure 2';
M2.metadata.info = 'UTCs 1-30, fields gamma & delta';

m2_inds = 1:30;
M2.Data.UTC = dstr;
M2.Data.UTC.Values = utcs(m2_inds);
M2.Data.gamma = dstr;
M2.Data.gamma.Values = gamma_vals(m2_inds);
M2.Data.delta = dstr;
M2.Data.delta.Values = delta_vals(m2_inds);

M3.metadata = mdata;
M3.metadata.file = 'tfile3.ict';
M3.metadata.description = 'Testing merge structure 3';
M3.metadata.info = 'UTCs 31-60, fields alpha & beta';

m3_inds = 31:60;
M3.Data.UTC = dstr;
M3.Data.UTC.Values = utcs(m3_inds);
M3.Data.alpha = dstr;
M3.Data.alpha.Values = alpha_vals(m3_inds);
M3.Data.beta = dstr;
M3.Data.beta.Values = beta_vals(m3_inds);

M4.metadata = mdata;
M4.metadata.file = 'tfile4.ict';
M4.metadata.description = 'Testing merge structure 4';
M4.metadata.info = 'UTCs 31-60, fields gamma & delta';

m4_inds = 31:60;
M4.Data.UTC = dstr;
M4.Data.UTC.Values = utcs(m4_inds);
M4.Data.gamma = dstr;
M4.Data.gamma.Values = gamma_vals(m4_inds);
M4.Data.delta = dstr;
M4.Data.delta.Values = delta_vals(m4_inds);

M1M2.metadata = mdata;
M1M2.metadata.file = {M1.metadata.file, M2.metadata.file};
M1M2.metadata.description = {M1.metadata.description, M2.metadata.description};
M1M2.metadata.info = {M1.metadata.info, M2.metadata.info};
M1M2.metadata.upper_lod_value = {M1.metadata.upper_lod_value, M2.metadata.upper_lod_value};
M1M2.metadata.lower_lod_value = {M1.metadata.lower_lod_value, M2.metadata.lower_lod_value};

M1M2.Data.UTC = dstr;
M1M2.Data.alpha = dstr;
M1M2.Data.beta = dstr;
M1M2.Data.gamma = dstr;
M1M2.Data.delta = dstr;
M1M2.Data.UTC.Values = utcs(sort(unique([m1_inds, m2_inds])));
M1M2.Data.alpha.Values = alpha_vals(m1_inds);
M1M2.Data.beta.Values = beta_vals(m1_inds);
M1M2.Data.gamma.Values = gamma_vals(m2_inds);
M1M2.Data.delta.Values = delta_vals(m2_inds);

M1M3.metadata = mdata;
M1M3.metadata.file = {M1.metadata.file, M3.metadata.file};
M1M3.metadata.description = {M1.metadata.description, M3.metadata.description};
M1M3.metadata.info = {M1.metadata.info, M3.metadata.info};
M1M3.metadata.upper_lod_value = {M1.metadata.upper_lod_value, M3.metadata.upper_lod_value};
M1M3.metadata.lower_lod_value = {M1.metadata.lower_lod_value, M3.metadata.lower_lod_value};

M1M3.Data.UTC = dstr;
M1M3.Data.alpha = dstr;
M1M3.Data.beta = dstr;
M1M3.Data.UTC.Values = utcs(sort(unique([m1_inds, m3_inds])));
M1M3.Data.alpha.Values(m1_inds) = alpha_vals(m1_inds);
M1M3.Data.beta.Values(m1_inds) = beta_vals(m1_inds);
M1M3.Data.alpha.Values(m3_inds) = alpha_vals(m3_inds);
M1M3.Data.beta.Values(m3_inds) = beta_vals(m3_inds);

M1M4.metadata = mdata;
M1M4.metadata.file = {M1.metadata.file, M4.metadata.file};
M1M4.metadata.description = {M1.metadata.description, M4.metadata.description};
M1M4.metadata.info = {M1.metadata.info, M4.metadata.info};
M1M4.metadata.upper_lod_value = {M1.metadata.upper_lod_value, M4.metadata.upper_lod_value};
M1M4.metadata.lower_lod_value = {M1.metadata.lower_lod_value, M4.metadata.lower_lod_value};

M1M4.Data.UTC = dstr;
M1M4.Data.alpha = dstr;
M1M4.Data.beta = dstr;
M1M4.Data.gamma = dstr;
M1M4.Data.delta = dstr;
M1M4.Data.UTC.Values = utcs(sort(unique([m1_inds, m4_inds])));
M1M4.Data.alpha.Values(m1_inds) = alpha_vals(m1_inds);
M1M4.Data.beta.Values(m1_inds) = beta_vals(m1_inds);
M1M4.Data.gamma.Values(m4_inds) = gamma_vals(m4_inds);
M1M4.Data.delta.Values(m4_inds) = delta_vals(m4_inds);


TotalMerge.metadata = mdata;
TotalMerge.metadata.file = {M1.metadata.file, M2.metadata.file, M3.metadata.file, M4.metadata.file};
TotalMerge.metadata.description = {M1.metadata.description, M2.metadata.description, M3.metadata.description, M4.metadata.description};
TotalMerge.metadata.info = {M1.metadata.info, M2.metadata.info, M3.metadata.info, M4.metadata.info};
TotalMerge.metadata.upper_lod_value = {M1.metadata.upper_lod_value, M2.metadata.upper_lod_value, M3.metadata.upper_lod_value, M4.metadata.upper_lod_value};
TotalMerge.metadata.lower_lod_value = {M1.metadata.lower_lod_value, M2.metadata.lower_lod_value, M3.metadata.lower_lod_value, M4.metadata.lower_lod_value};

TotalMerge.Data.UTC = dstr;
TotalMerge.Data.alpha = dstr;
TotalMerge.Data.beta = dstr;
TotalMerge.Data.gamma = dstr;
TotalMerge.Data.delta = dstr;
TotalMerge.Data.UTC.Values = utcs;
TotalMerge.Data.alpha.Values = alpha_vals;
TotalMerge.Data.beta.Values = beta_vals;
TotalMerge.Data.gamma.Values = gamma_vals;
TotalMerge.Data.delta.Values = delta_vals;

end