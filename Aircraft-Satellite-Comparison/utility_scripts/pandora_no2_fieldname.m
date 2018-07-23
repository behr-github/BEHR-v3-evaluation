function no2_field = pandora_no2_fieldname(Merge)
%PANDORA_NO2_FIELDNAME Return the correct NO2 field for a given Merge struct
%
%   NO2_FIELD = PANDORA_NO2_FIELDNAME( MERGE ) Finds the NO2 field 

E = JLLErrors;

if isfield(Merge.Data, 'NO2')
    no2_field = 'NO2';
elseif isfield(Merge.Data, 'NO2Molcm2')
    no2_field = 'NO2Molcm2';
else
    E.callError('unknown_no2_field', 'Cannot identify NO2 field in %s', Merge.metadata.file);
end

if strcmpi(Merge.Data.(no2_field).Unit, {'DU','Dobson Unit','Dobson Units'})
    E.callError('wrong_no2_field', 'The %s field in %s is in Dobson Units, not molecules cm^-2', no2_field, Merge.Data.(no2_field).Unit)
elseif ~any(strcmpi(Merge.Data.(no2_field).Unit, {'mol/cm2', 'molecules/cm2'}))
    warning('pandora_no2_field:unit_not_recognized','NO2 units "%s" not recognized as being molecules cm^-2', Merge.Data.(no2_field).Unit)
end

end

