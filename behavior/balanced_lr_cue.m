function [balanced_left, balanced_right] = balanced_lr_cue(trials,locs);
% locs are either 11 or 12
% trials is logical

trials_leftcued = (trials.*(locs==11))==1;
trials_rightcued = (trials.*(locs==12))==1;
left_cued = find(trials_leftcued);
right_cued = find(trials_rightcued);
n_left = length(left_cued);
n_right = length(right_cued);
n_equal = min([n_left,n_right]);
left_cued = randsample(left_cued,n_equal);
right_cued = randsample(right_cued,n_equal);
balanced_left = zeros(size(trials_leftcued));
balanced_right= zeros(size(trials_leftcued));
balanced_left(left_cued)=1;
balanced_right(right_cued)=1;

end