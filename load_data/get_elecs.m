function [elecs,lfp_elec] = get_elecs(session)

if strcmp(session,'pmc_gaba1')
    elecs = 1:4;
    lfp_elec = 4;
elseif strcmp(session,'pmc_gaba2')
    elecs = 1:4;
    lfp_elec = 4;
elseif strcmp(session,'pmc_sal1')
    elecs = 1:4;
    lfp_elec = 4;
elseif strcmp(session,'pmc_sal2')
    elecs = 1:4;
    lfp_elec = 4;
elseif strcmp(session,'pmc_sal3')
    elecs = 1:4;
    lfp_elec = 4;
elseif strcmp(session,'pmc_gaba3')
    elecs = 5:8;
    lfp_elec = 8;
elseif strcmp(session,'pmc_gaba4')
    elecs = 5:8;
    lfp_elec = 8;
elseif strcmp(session,'put_gaba')
    elecs = 1:4;
    lfp_elec = 2;
elseif strcmp(session,'put_noic1')
    elecs = 1:4;
    lfp_elec = 3;
elseif strcmp(session,'put_noic2')
    elecs = 1:4;
    lfp_elec = 3;
elseif strcmp(session,'put_noic3')
    elecs = 1:4;
    lfp_elec = 3;
elseif strcmp(session,'pmc_acsf1')
    elecs = 29:32;
    lfp_elec = 29;
elseif strcmp(session,'pmc_acsf2')
    elecs = [];
    lfp_elec = [];
elseif strcmp(session,'pmc_noic1')
    elecs = 1:2:7;
    lfp_elec = 1;
elseif strcmp(session,'pmc_noic2')
    elecs = 1:2:7;
    lfp_elec = 1;
elseif strcmp(session,'pmc_noic3')
    elecs = 1:2:7;
    lfp_elec = 1;
elseif strcmp(session,'pmc_noic4')
    elecs = 1:2:7;
    lfp_elec = 1;
end

end
