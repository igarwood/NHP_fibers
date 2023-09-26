%% Plot fiber electrode impedance before and after autoclave

meta_data = setup;
folder = [meta_data.device_folder,'impedance/'];

devs = [1,2,3,4,5];
elecs = 1:4;
dev_prefix = 'dev';
elec_prefix = '_elec';
elec_desc = '';%_postautoclave

z_test_post = [];
test_freq = 1000;

z_all = zeros(31,length(devs)*4);
for i = 1:length(devs)
    dev = devs(i);
    elecs_i = elecs;
    for j = 1:length(elecs_i)
        elec = elecs_i(j);
        z_mat = readmatrix([folder,dev_prefix,num2str(dev),...
            elec_desc,elec_prefix,num2str(elec)]);
        [~,z_test_ind] = min(abs(z_mat(:,1)-test_freq));
        z_vec = sqrt(sum(z_mat(:,2:3).^2,2));
        z_all(:,(i-1)*4+j) = z_vec;
        z_test_post = [z_test_post, sqrt(sum(z_mat(z_test_ind,2:3).^2))];
    end
end

z = mean(z_test_post)/1000;
z_sd = std(z_test_post)/1000;

z_all = z_all/1000;


mean_z = mean(z_all,2);
CI_z = zeros(2,length(mean_z));
CI_z(1,:) = mean_z-std(z_all,[],2)/sqrt(size(z_all,2));
CI_z(2,:) = mean_z+std(z_all,[],2)/sqrt(size(z_all,2));
se = std(z_all,[],2)/sqrt(size(z_all,2));
freqs = z_mat(:,1)';

defaults = [
        0.4940    0.1840    0.5560;...
        0.0000    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.6350    0.0780    0.1840];

figure
x =[freqs, fliplr(freqs)];
y =[CI_z(1,:), fliplr(CI_z(2,:))];

plot(freqs,mean_z,'.','Color',defaults(1,:),'HandleVisibility','off');
hold on
fill(x,y,defaults(2,:),'FaceAlpha',0.3,'EdgeColor',...
    defaults(2,:));


devs = [5];
elecs = 1:4;
elec_desc = '_postAC';

z_test_post = [];
test_freq = 1000;

z_all = zeros(31,length(devs)*4);
for i = 1:length(devs)
    dev = devs(i);
    elecs_i = elecs;
    for j = 1:length(elecs_i)
        elec = elecs_i(j);
        z_mat = readmatrix([folder,dev_prefix,num2str(dev),...
            elec_desc,elec_prefix,num2str(elec)]);
        [~,z_test_ind] = min(abs(z_mat(:,1)-test_freq));
        z_vec = sqrt(sum(z_mat(:,2:3).^2,2));
        z_all(:,(i-1)*4+j) = z_vec;
        z_test_post = [z_test_post, sqrt(sum(z_mat(z_test_ind,2:3).^2))];
    end
end

z = mean(z_test_post);
z_sd = std(z_test_post);

z_all = z_all/1000;

mean_z = mean(z_all,2);
CI_z = zeros(2,length(mean_z));
CI_z(1,:) = mean_z-std(z_all,[],2)/sqrt(size(z_all,2));
CI_z(2,:) = mean_z+std(z_all,[],2)/sqrt(size(z_all,2));
se = std(z_all,[],2)/sqrt(size(z_all,2));


y =[CI_z(1,:), fliplr(CI_z(2,:))];

plot(freqs,mean_z,'.','Color',defaults(1,:),'HandleVisibility','off');
hold on
fill(x,y,defaults(3,:),'FaceAlpha',0.3,'EdgeColor',...
    defaults(3,:));

set(gca,'XScale','log')
set(gcf,'renderer','painters')
xlabel('Frequency (Hz)')
ylabel('Impedance (kOhm)')

legend('Before autoclave','After autoclave')
