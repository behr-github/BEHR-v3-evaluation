function [ temperature_prof, pres_profile ] = campaign_temperature_profile( campaign_name )
%campaign_temperature_profile Returns the median temp profile for a given campaign.
%   Over the course of a typical field campaign, it is unlikely that the
%   free troposphere temperature will vary that significantly, so these
%   profiles were derived from temperature data for entire campaign (see
%   make_composite_profile.m in NO2 Profiles/One-off scripts).

E = JLLErrors;

if ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'md'))
    temperature_prof = [300.4550 299.9400 298.9150 298.9100 299.1000 299.0700 296.6100 295.0600 293.5300 290.4100 290.9900 287.9900 286.7000 285.6500 284.5300 281.0800 280.9000 276.4500 272.8400 268.0562 263.6088 258.6371 253.0005 243.5814 229.3784];
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'ca'))
    temperature_prof = [283.5403 284.4901 283.9849 284.1614 294.1300 296.7400 287.3189 292.4100 288.3163 287.3151 288.9400 286.2400 283.9400 282.3737 279.8776 280.3500 280.7500 276.2100 272.7200 251.4486 251.5334 251.6282 251.7357 251.9153 252.1861];
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'tx'))
    temperature_prof = [296.5900 299.4791 288.6492 287.3828 297.9100 297.4400 293.9000 294.5200 291.6400 288.9400 289.9278 287.5900 285.8800 284.1500 280.3100 281.1402 280.0400 276.4336 272.7800 251.4486 251.5334 251.6282 251.7357 251.9153 252.1861];
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'co'))
    temperature_prof = [300.3389 300.1105 299.8799 299.6469 299.2929 298.9335 298.5684 298.0725 297.4374 296.7843 296.1123 294.7073 293.9821 291.2372 289.3689 285.8652 282.4206 277.2497 271.9121 268.6394 260.4381 251.2698 240.8757 223.5061 197.3150];
elseif ~isempty(regexpi(campaign_name,'seac4rs')) || ~isempty(regexpi(campaign_name,'seacers'));
    temperature_prof = [300.6500 300.2505 297.6500 288.3628 298.1900 297.6500 296.1088 294.9000 293.0738 290.9000 290.3594 288.6500 286.6884 285.1200 280.7724 281.5697 279.3009 275.6400 272.0000 265.1500 262.6500 254.1500 248.9000 237.9000 224.9000];
elseif ~isempty(regexpi(campaign_name,'intex')) && ~isempty(regexpi(campaign_name,'b'))
    temperature_prof = [300.2849 287.8000 289.8130 289.1352 297.9000 297.4000 295.9000 294.7948 292.7500 290.8000 290.1500 288.4824 286.3800 284.9700 280.9400 281.4000 279.0000 274.9000 270.9000 263.4000 261.0000 252.2000 243.5000 232.4000 221.1500];
elseif ~isempty(regexpi(campaign_name, 'arctas')) && ~isempty(regexpi(campaign_name, 'carb'))
    temperature_prof = [296.426354295365 294.05 291.65 293.65 300.65 301.3999939 299.95 297.65 298.1499939 294.35 293.6499939 293.6499939 294.3999939 290.15 288.35 282.6499939 279.95 273.8999939 270.6499939 263.6499939 256.6499939 249.3999939 244.15 232.8999939 223.6499939];
elseif ~isempty(regexpi(campaign_name, 'dc3'))
    temperature_prof = [296.23978412048 296.875244511785 297.517027995884 298.165261673296 299.1499939 300.1499939 294.3999939 294.1499939 296.8999939 293.3999939 293.1499939 293.8999939 294.6499939 291.1499939 288.3999939 282.3999939 279.8999939 273.8999939 270.6499939 263.3999939 256.6499939 249.3999939 244.1499939 231.6499939 223.6499939]
elseif ~isempty(regexpi(campaign_name,'soas'))
    temperature_prof = [301.3045 300.359 302.032 299.124 296.399 298.121 295.871 296.134 294.472 290.12 289.809 289.256 286.562 285.327 282.8805 282.16 278.988 274.693 272.115 267.025 264.347 261.3532 257.9592 252.2875 243.7352];
else
    error(E.badinput('Could not parse the given campaign name'));
end

pres_profile = bin_omisp_pressure(1,1);
if temperature_prof(end) > 230
    warning('Temperature profile does not reach the expected range of temperatures for the tropopause. Use with caution.');
end
    
end

