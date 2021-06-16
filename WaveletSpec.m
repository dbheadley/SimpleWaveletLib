function spec = WaveletSpec(data, tVal, freqs, varargin)
  % NOTE: MULTIPLE CHANNEL SUPPORT IS BEING DEPRECATED, DO NOT USE
  % data is a vector of samples across time (i.e. time series), or an array
  % with samples along the column axis and channels/epochs along the row axis.
  % tVal can either be a scalar, in which case it is the the time between 
  % samples, or a vector with the same length as data, in which case it is
  % the corresponding time for each sample. 
  % freqs are the frequencies you want to calculate your spectrum for in Hz
  
  % Optional Arguments
  % 'SAVE' - followed by the file path to save the spectral output to. A
  % default file suffix of .wls is attached if not already provided.
  % For the output, instead of the normal spectrogram structure being 
  % returned, only the filepath where it can be found is. 
  % 'WAVELETSPEC' - followed by a 2x1 cell array of wavelet parameters. The
  % first entry is a string specifying the wavelet type, by default this is
  % 'morlet'. The second is a list of parameters specific for that wavelet
  % type. For the morlet wavelet the only parameter is sigma.
  % 'CHANNELNAME' - followed by a character array or cell array of strings,
  % the name of the channel on which the spectrogram is calculated. 
  % By default channels are named 'Ch1' to 'ChN', where N is the number of 
  % rows in data.
  % 'SUBSAMPLE' - followed by an integer, downsamples the returned wavelet
  % spectrogram by a user specified factor e.g. setting factor to 10 will
  % return every 10th time point.
  % 'OVERWRITE' - if specified when SAVE is used, then allows for previously created wls
  % files to be overwritten. Otherwise, the function just returns the already 
  % existing spec.
  % 'NORMALIZE' - if specified, then will use the appropriate normalization
  % for a given wavelet type.
  % 'PAD' - if specified, then the signal is reflected at either end to
  % mitigate edge effects.
  
  % for compatibility with the multichannel case, make data a row vector 
  % if supplied as a column vector
  if size(data,1)==numel(data)
    data = data';
  end
  
  numChans = size(data,1);
  numSamps = size(data,2);
  
  if any(strcmp(varargin,'SAVE'))
    fileName = varargin{find(strcmp(varargin,'SAVE'))+1};
    if ~strcmp(fileName((end-3):end),'.wls')
      fileName = [fileName '.wls'];
    end
    saveYes = true;
    % if file already exists and overwrite is not specified then exit with
    % file name
    if any(strcmp(varargin, 'OVERWRITE'))
      overWriteYes = true;
    else
      overWriteYes = false;
    end
  else
    saveYes = false;
  end
  
  if any(strcmp(varargin,'CHANNELNAME'))
    channelNames = varargin{find(strcmp(varargin,'CHANNELNAME'))+1};
    if ischar(channelNames)
      channelNames = {channelNames};
    end
    if length(channelNames) ~= numChans
      error('The wrong number of channel names were specified');
    end
  else
    for j = 1:numChans
      channelNames{j} = ['Ch' num2str(j)];
    end
  end
  
  if any(strcmp(varargin, 'WAVELETSPEC'))
    wType = varargin{find(strcmp(varargin,'WAVELETSPEC'))+1}{1};
    wParams = varargin{find(strcmp(varargin,'WAVELETSPEC'))+1}{2};
  else
    wType = 'morlet';
    wParams = 7;
  end
  
  if any(strcmp(varargin, 'NORMALIZE'))
    normYes = true;
  else
    normYes = false;
  end
  
  
  
  if any(strcmp(varargin, 'SUBSAMPLE'))
    subFactor = varargin{find(strcmp(varargin,'SUBSAMPLE'))+1};
    retInds = 1:subFactor:numSamps;
  else
    retInds = 1:numSamps;
  end
  
  if any(strcmp(varargin, 'PAD'))
    padYes = true;
  else
    padYes = false;
  end
  
  if length(tVal) == 1
    spec.Times = (0:(numSamps-1))*tVal;
    spec.Times = spec.Times(retInds);
    tStep = tVal;
  elseif length(data) == length(tVal)
    spec.Times = tVal(retInds);
    tStep = median(diff(tVal));
  else
    error('Time improperly specified');
  end
  
  if padYes
    retInds = retInds+size(data,2);
    data = [data data data];
  end
  
  if saveYes
    sampleTimes = spec.Times;
    Spectrum = 1i*ones(1,numel(retInds),numChans);%numel(freqs),numel(spec.Times));
    waveletParams = {wType wParams};
    
    if isempty(dir(fileName)) || overWriteYes
      save(fileName,'sampleTimes', 'freqs', 'Spectrum', 'channelNames', ...
        'waveletParams', '-v7.3');
      clear sampleTimes Spectrum waveletParams spec;
      mfObj = matfile(fileName, 'Writable',true);
      spec = fileName;
    elseif ~overWriteYes
      spec = fileName;
      return;
    end
  else
    spec.Freqs = freqs;
    spec.Spectrum = 1i*ones(numel(spec.Freqs),numel(retInds),numChans);
    spec.WaveletParams = {wType wParams};
    spec.ChannelNames = channelNames;
  end
  
  
  for j = 1:length(freqs)
    currFreq = freqs(j);
    
    switch wType
      case 'morlet'
        sigma = wParams(1)/(2*pi*currFreq);
        amp = 1/(sigma*sqrt(2*pi));
        tSupp = (-wParams(1)/(2*currFreq)):tStep:(wParams(1)/(2*currFreq));
        cmwl = amp*exp(-(tSupp.^2)./(2*(sigma^2))).*exp(2*1i*pi*currFreq*tSupp);
        if normYes
          cmwl = sqrt((2*pi*(1/currFreq))/tStep)*cmwl;
        end
        w1 = real(cmwl);
        w2 = imag(cmwl);
      otherwise
        error('No valid wavelet');
    end
    
    if numChans == 1
      if saveYes
        temp = conv(data(1,:),w1,'same')+1i*conv(data(1,:),w2,'same');
        mfObj.Spectrum(j,:) = temp(retInds);
      else
        temp = conv(data(1,:),w1,'same')+1i*conv(data(1,:),w2,'same');
        spec.Spectrum(j,:) = temp(retInds);
      end
    else
      for k = 1:numChans
        if saveYes
          temp = conv(data(k,:),w1,'same')+1i*conv(data(k,:),w2,'same');
          mfObj.Spectrum(j,:,k) = temp(retInds);
        else
          temp = conv(data(k,:),w1,'same')+1i*conv(data(k,:),w2,'same');
          spec.Spectrum(j,:,k) = temp(retInds);
        end
      end
    end
  end
    