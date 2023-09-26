function modulation_trials = get_modulation_trials(trial_times,...
    modulation_start,modulation_end);


modulation_start_trial = zeros(size(modulation_start));
modulation_end_trial = zeros(size(modulation_end));
try
    for n = 1:length(modulation_start)
       modulation_start_trial(n) = ...
           find(trial_times(2,:) > modulation_start(n),1);
       modulation_end_trial(n) = ...
           find(trial_times(1,:) >= modulation_end(n),1)-1;
    end

    modulation_trials = [modulation_start_trial,modulation_end_trial];
catch
    modulation_trials = [];
end

end