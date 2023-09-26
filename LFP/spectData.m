function [t, stimes, sfreqs, time_series, spectral_data] = ...
    spectData(time_series,varargin)
% Multitaper spectral analysis of time_series data

defaults.moving_win = [30 10];
defaults.tapers = [10 19];

if any(strcmp(varargin,'time'))
    t = varargin{find(strcmp(varargin,'time'))+1};
else
    if any(strcmp(varargin,'fs'))
        t = (0:(length(time_series)-1))*(1/fs);
    else
        error('spectData input must include time or fs')       
    end  
end

if any(strcmp(varargin,'fs'))
    fs = varargin{find(strcmp(varargin,'fs'))+1};
else
    fs = 1/(t(2)-t(1))
end

if any(strcmp(varargin,'moving_win'))
    moving_win = varargin{find(strcmp(varargin,'moving_win'))+1};
else
    moving_win = defaults.moving_win;
end

if any(strcmp(varargin,'tapers'))
    tapers = varargin{find(strcmp(varargin,'tapers'))+1};
else
    tapers = defaults.tapers;
end


time_series = time_series-mean(time_series);
params.Fs=fs;
params.tapers=tapers;

% Estimate the spectrogram using mutitaper spectral estimation
[spectral_data, stimes, sfreqs]=mtspecgramc(time_series, moving_win, params);

end