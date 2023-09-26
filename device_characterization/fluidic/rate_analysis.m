%% Extract and plot fiber fluidic flow rate data

meta_data = setup;
folder = [meta_data.device_folder,'fluidic/'];

filebase = {'10_nlm_trial_',...
    '20_nlm_trial_',...
    '30_nlm_trial_',...
    '40_nlm_trial_',...
    '50_nlm_trial_',...
    '60_nlm_trial_',...
    '70_nlm_trial_',...
    '80_nlm_trial_',...
    '90_nlm_trial_',...
    '100_nlm_trial_'};

ntrials = 3;
nrates = size(filebase,2);
trials = [1,2,3];

%% Extract data
trial_time = cell(nrates,ntrials);
trial_volume = cell(nrates,ntrials);

for k = 1:nrates
    file_k = filebase{k};
    for n = 1:ntrials
        trial_n = trials(n);
        load([folder,file_k,num2str(trial_n)]);
        volume = volume*1e-6;
        trial_time{k,n} = time;
        trial_volume{k,n} = volume;
    end
end

%% Align and estimate rate

start_time_trial = zeros(1,ntrials);
trial_time_aligned = cell(nrates,ntrials);
trial_volume_aligned = cell(nrates,ntrials);
trials_time = cell(nrates,1);
trials_volume = cell(nrates,1);

% Remove timepoints beyond the end of the infusion
max_time = [1400,650,400,300,350,320,240,250,220,220];
rate = zeros(1,10);
rate_se = zeros(1,10);

for k = 1:nrates
    for n = 1:ntrials
        start_time_trial(n) = trial_time{k,n}(1);
    end
    start_time = max(start_time_trial);
    for n = 1:ntrials
        ind = ((trial_time{k,n}>start_time).*...
            (trial_time{k,n}<max_time(k))) == 1;
        trial_time_aligned{k,n} = trial_time{k,n}(ind);
        trial_volume_aligned{k,n} = trial_volume{k,n}(ind);
        trial_volume_aligned{k,n} = trial_volume_aligned{k,n} - ...
            trial_volume_aligned{k,n}(1);
    end
    
    trials_time{k} = [trial_time_aligned{k,1},...
        trial_time_aligned{k,2},...
        trial_time_aligned{k,3}];
    trials_volume{k} = [trial_volume_aligned{k,1};...
        trial_volume_aligned{k,2};...
        trial_volume_aligned{k,3}];
    
    mdl = fitlm(trials_time{k}/60,trials_volume{k});
    coeff = table2array(mdl.Coefficients);
    rate(k) = coeff(2,1);
    rate_se(k) = coeff(2,2);
end
%% Plot observed rate vs set rate
rates = 10:10:100;
figure
plotConfInterval(rates,rate,rate-rate_se*1.96,rate+rate_se*1.96,[],20)
hold on
rates = 0:10:110;
plot(rates,rates,'k--')
xlabel('Set infusion rate (nl/min');
ylabel('Measured infusion rate (nl/min)');
xlim([0,110])
ylim([0,110])

%% Stats

sigma_hat_sq = mdl.NumObservations^(-1)*sum((mdl.Variables{:,2}-...
    mdl.Coefficients{1,1}-mdl.Coefficients{2,1}*mdl.Variables{:,1}).^2);

% Equal to mdl.MSE:
s_hat_sq = mdl.NumObservations*(mdl.NumObservations-2)^(-1)*sigma_hat_sq;

% Equal to mdl.Coefficients SE
B_error = sqrt(s_hat_sq/sum((mdl.Variables{:,1}-...
    mean(mdl.Variables{:,1})).^2));


%% Plot example measured vs set rate data

mdl = fitlm(trials_time{5}/60,trials_volume{5});

x = linspace(min(trials_time{5}/60),max(trials_time{5}/60));
y_error = 1.96*sqrt(1/mdl.NumObservations + (x-...
    mean(mdl.Variables{:,1})).^2)/sum((x-...
    mean(mdl.Variables{:,1})).^2);
y_est = predict(mdl,x');
y_set = x*50+mdl.Coefficients{1,1};
figure
plot(trials_time{5},trials_volume{5},'.')
hold on
plot(x*60,y_set,'k--')
xlim([min(trials_time{5}),max(trials_time{5})])
xlabel('Time (sec)')
ylabel('Volume (nl)')
