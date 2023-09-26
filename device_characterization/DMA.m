%% Plot dynamic materials analysis results

meta_data = setup;
folder = [meta_data.device_folder,'DMA/'];

filenames = {'fiber_y.csv','steel.csv'};
y_labels = [1,2,3];
n_samples = [3,1];

%% Extract data
samples_remove = cell(1,2);
samples_remove{1} = [];
load('color_defaults');

labels = {'fiber','steel'};
N = size(n_samples,2);
J = max(n_samples);
DMA_data= cell(J,N);
all_DMA_data = cell(1,N);


for n = 1:N
    for j = 1:n_samples(n)
        if ~ismember(j,samples_remove{n})
            ind_sample = strfind(filenames{n},'y');
            if ~isempty(ind_sample)
                filename = [filenames{n}(1:ind_sample-1),...
                    num2str(y_labels(j)),filenames{n}(ind_sample+1:end)];
            else
                filename = filenames{n};
            end
            
            T = readtable([folder,filename]);
            [freq_row,freq_col] = find(strcmp(['Frequency'],T{:,:}));
            [stiffness_row,stiffness_col] = ...
                find(strcmp(['Stiffness'],T{:,:}));
            DMA_data_str = T{[(freq_row+2):end],...
                [freq_col,stiffness_col]};
            DMA_data{j,n} = zeros(size(DMA_data_str));
            for a = 1:size(DMA_data_str,1)
                for b = 1:size(DMA_data_str,2)
                    DMA_data{j,n}(a,b) = str2num(DMA_data_str{a,b});
                end
            end
            all_DMA_data{n} = [all_DMA_data{n},DMA_data{j,n}(:,2)];

        end
    end
    
end


%% Plot results
figure
DMA_fiber = mean(all_DMA_data{1}(2:end,:),2);
DMA_fiber_std = std(all_DMA_data{1}(2:end,:),[],2);
DMA_fiber_se = DMA_fiber_std/sqrt(size(all_DMA_data{1}(2:end,:),2));
DMA_fiber_ci = [DMA_fiber+DMA_fiber_se, DMA_fiber-DMA_fiber_se];

freq = DMA_data{1,1}(2:end,1);
plotConfInterval(freq',DMA_fiber',DMA_fiber_ci(:,1)',...
    DMA_fiber_ci(:,2)');
hold on

plot(freq,all_DMA_data{2}(2:end),'.','color',defaults(2,:))


set(gca, 'XScale', 'log')
xlim([10^-2,10])
xlabel('Frequency (Hz)')
ylabel('Stiffness (N/m)')
[~,f_end] = min(abs(freq-10)); 
ylim([0,200])
DMA_fiber_all = reshape(all_DMA_data{1}(2:f_end,:),1,[]);
DMA_fiber_all_mean = mean(DMA_fiber_all);
DMA_fiber_all_std = std(DMA_fiber_all);
DMA_fiber_all_se = DMA_fiber_all_std/sqrt(length(DMA_fiber_all));

DMA_steel_all = reshape(all_DMA_data{2}(2:f_end,:),1,[]);
DMA_steel_all_mean = mean(DMA_steel_all);
DMA_steel_all_std = std(DMA_steel_all);
DMA_steel_all_se = DMA_steel_all_std/sqrt(length(DMA_steel_all));


% t-test 
DMA_diff = all_DMA_data{1}(2:f_end,:)-all_DMA_data{2}(2:f_end,:);
[h,p] = ttest(reshape(DMA_diff,1,[]))
