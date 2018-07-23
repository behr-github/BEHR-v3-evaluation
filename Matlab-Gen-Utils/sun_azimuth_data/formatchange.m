function [degree_part,minute_part,second_part] = formatchange(decimalvalue)

degree_part = floor(decimalvalue);
fraction_part = decimalvalue - floor(decimalvalue);
minute_part = floor(fraction_part * 60);
second_part = (fraction_part * 3600) - (minute_part * 60);