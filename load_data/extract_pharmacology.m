function [properties] = extract_pharmacology(file)

load(file,'infusion');
if isempty(infusion)
    properties = [];
else
    properties.infusion_start = infusion.start;
    properties.infusion_end = infusion.end;
    properties.infusion_rate = infusion.rate;
    properties.drug = infusion.drug;
end

end