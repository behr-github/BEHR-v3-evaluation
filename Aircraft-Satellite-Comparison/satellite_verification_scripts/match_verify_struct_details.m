function details_struct = match_verify_struct_details(comp_struct, details_struct)
%MATCH_VERIFY_STRUCT_DETAILS Cut down the details field to match the others
%   VSTRUCT = MATCH_VERIFY_STRUCT_DETAILS( COMP_STRUCT, DETAILS_STRUCT )
%   The first structure returned by run_insitu_verification() has the main
%   fields that are vectors of VCDs and other data for valid profile
%   comparisons, and the second is a structure of all profiles considered.
%   This is useful if you need to check why a particular profile was
%   omitted, but is inconvenient when trying to match particular VCD
%   comparisons against detailed profile information. This function takes
%   in one of those pairs of structures and modifies the details one so
%   that indexing it returns the details that match any of the fields in
%   the main structure at the same index.

E = JLLErrors;

details_prof_ids = {details_struct.profile_id};
details_perm_vec = nan(size(comp_struct.profile_ids));
for i_prof = 1:numel(comp_struct.profile_ids)
    this_prof_id = comp_struct.profile_ids{i_prof};
    xx_detail = find(cellfun(@(x) isequal(x, this_prof_id), details_prof_ids));
    if isempty(xx_detail)
        E.callError('details_not_found', 'Could not find details for profile id %s', num2str(this_prof_id));
    elseif numel(xx_detail) > 1
        E.callError('multiple_details_found', 'Found multiple details (at %s) for profile id %s', num2str(xx_detail), num2str(this_prof_id));
    else
        details_perm_vec(i_prof) = xx_detail;
    end
end

details_struct = details_struct(details_perm_vec);

end

