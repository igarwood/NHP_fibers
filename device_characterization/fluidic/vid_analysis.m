%% Extract volume infused through capillary from raw video data

meta_data = setup;
folder = [meta_data.device_folder,'fluidic/raw_videos/'];
example_file = '50nlm_oil_trial1.mp4';
filename = [folder,example_file];

trial = '50_nlm_trial_1';
save_data = 0; % set to 1 to save volume and time data

start_time = 10;
% This might not work on linux. Try converting to .avi or updating
% gstreamer:
vid = VideoReader(filename); 

% Video meta deta
frame = read(vid,1);
frame_gray = rgb2gray(frame);
frame_double = im2double(frame_gray);
figure
imagesc(frame_double) 
width = vid.Width;
height = vid.Height;


%% Extract first video frame

slice = frame_double(width/2,:);
diff_slice = diff(slice);
[~,left_pix] = min(diff_slice);
[~,right_pix] = max(diff_slice);

figure
imagesc(frame_double)
hold on
plot([left_pix, left_pix],[0, height],'w','linewidth',3)
plot([right_pix, right_pix],[0, height],'w','linewidth',3)


%% Capillary dimensions in pixels and um
cap_od_pix = right_pix-left_pix;
cap_od = 375;
cap_id = 100;
length_scale = cap_od/cap_id;


xarea = (left_pix+50):(right_pix-50);
xarea_width = xarea(end)-xarea(1);

%% Extract each frame and compute volume over time    
thresh = 0.8;

nframes = vid.NumFrames;
boundary_start = zeros(nframes,1);
for ind = 1:nframes
    frame = read(vid,ind);
    frame_gray = rgb2gray(frame);
    frame_double = im2double(frame_gray);
   
    
    area = frame_double(:,xarea);
    area_thresh = area<thresh;
    boundary = sum(area_thresh,2);
    boundary = boundary>(xarea_width*0.9);
    start_0 = find(diff(boundary)==-1);
    if length(start_0) == 1
        boundary_start(ind) = start_0;
    elseif length(start_0) > 1
        boundary_start(ind) = max(start_0);
    end
end

diff_boundary = diff(boundary_start);
switch_ind = find(diff_boundary > 100);

time = [0:(nframes-1)]./vid.FrameRate-start_time;
time(boundary_start == 0) = [];
boundary_start(boundary_start == 0) = [];

distance = boundary_start*length_scale;
area = pi*(cap_id/2)^2;
volume = area*distance;
volume = -(volume-volume(1));

figure
plot(time,volume*1e-6,'.')
xlabel('time (s)')
ylabel('volume (ul)')
hold on

if save_data==1
    save(trial,'time','volume');
end

