classdef WaveletObj < dynamicprops
  % Acts as an interface to a wavelet spectrogram data file, which has the
  % extension .wls. You can only read from the file, you cannot modify it.
  
  properties (SetAccess = private, GetAccess = public)
    NumFreqs;
    NumSamps;
    NumChans; % DEPRECATED
    Duration;
    WaveletParams
    Freqs
    FileLocation
    DataTypes
    TimeStep
  end
  
  properties (SetAccess = private, GetAccess = private)
    mfObj
  end
  
  methods
    % Create a WaveletObj from an extant file
    function obj = WaveletObj(fileLocation, varargin)
      
      % confirm that file is valid
      if ~strcmp(fileLocation((end-3):end),'.wls')
        error('Wrong file type');
      end
      
      % if file does not exist, return an error
      if ~fileattrib(fileLocation)
        error('File does not exist');
      end
      
      % create matfile object and extract properties
      obj.mfObj = matfile(fileLocation,'Writable',false);
      varsPresent = whos('-file', fileLocation);
      % data variables are those that do not have the name 'channelNames',
      % 'freqs', 'sampleTimes', and 'waveletParams'
      
      % spectral analysis parameters
      obj.Freqs = obj.mfObj.freqs;
      obj.FileLocation = fileLocation;
      obj.WaveletParams = obj.mfObj.waveletParams;
      
      % spectral data
      dataVars = setdiff({varsPresent.name},{'channelNames' 'freqs' ...
        'waveletParams' 'sampleTimes'}); % channelNames is left here for legacy reasons
      obj.DataTypes = dataVars;
      
      
      % derived properties
      obj.TimeStep = median(diff(obj.mfObj.sampleTimes));
      obj.NumFreqs = numel(obj.Freqs);
      obj.NumSamps = numel(obj.mfObj.sampleTimes);
%       obj.NumChans = numel(obj.mfObj.channelNames);
      obj.Duration = obj.mfObj.sampleTimes(1,end);
    end
    
    % get spectral data
    function spec = LoadSpectrum(obj, varargin)
      spec = LoadWLSData(obj, 'Spectrum',varargin{:});
    end
    
    % get coherence data
    function spec = LoadCoherence(obj, varargin)
      spec = LoadWLSData(obj, 'Coherence',varargin{:});
    end
    
    % get phase data
    function spec = LoadPhase(obj, varargin)
      if any(strcmp(obj.DataTypes, 'Phase'))
        spec = LoadWLSData(obj, 'Phase',varargin{:});
      elseif any(strcmp(obj.DataTypes, 'Spectrum'))
        spec = LoadWLSData(obj, 'Spectrum', varargin{:});
        for j = 1:length(spec)
          spec(j).Phase = angle(spec(j).Spectrum);
        end
        spec = rmfield(spec,'Spectrum');
      else
        error('No data type that is compatible with phase');
      end
    end
    
    % get amplitude data
    function spec = LoadAmplitude(obj, varargin)
      if any(strcmp(obj.DataTypes, 'Amplitude'))
        spec = LoadWLSData(obj, 'Amplitude',varargin{:});
      elseif any(strcmp(obj.DataTypes, 'Spectrum'))
        spec = LoadWLSData(obj, 'Spectrum', varargin{:});
        for j = 1:length(spec)
          spec(j).Amplitude = abs(spec(j).Spectrum);
        end
        spec = rmfield(spec,'Spectrum');
      else
        error('No data type that is compatible with amplitude');
      end
    end
    
    % get sample time, emulate the behavior of the 'Times' field in the
    % spectrogram structure
    function tVals = Times(obj, varargin)
      if nargin == 1
        tVals = obj.mfObj.sampleTimes;
      elseif nargin == 2
        tVals = obj.mfObj.sampleTimes(1,varargin{1});
      elseif nargin == 3
        tVals = obj.mfObj.sampleTimes(varargin{1},varargin{2});
      else
        error('Time not specified correctly');
      end
    end
    
    % get sample index closest to time
    function sampleInd = GetTimeIndex(obj, tVal)
      if (tVal < obj.mfObj.sampleTimes(1,1)) || (tVal > obj.mfObj.sampleTimes(1,end))
        error('Time is out of range');
      end
      sampleInd = FindTimeIndex(obj,tVal);
    end
    
    % get sample index closest to frequency
    function sampleInd = GetFrequencyIndex(obj, fVal)
      if (fVal < obj.Freqs(1)) || (fVal > obj.Freqs(1,end))
        error('Frequency is out of range');
      end
      [~,sampleInd] = min(abs(fVal-obj.Freqs));
    end
    
  end
  
  methods (Access = private)
    
    tInd = FindTimeIndex(obj, tVal, varargin)
    
    function spec = LoadWLSData(obj, dType, varargin)
      % dType is a character array specifying the desired type of data to
      % return
      
      % test that specified data field is present
      if ~any(strcmp(obj.DataTypes,dType))
        error([dType ' not present']);
      end
      
      % if no inputs are given, then none of the spectrogram is returned
      if isempty(varargin)
        spec = struct('Freqs',[],'SampleTimes',[], dType, [], ...
          'WaveletParams', {obj.WaveletParams});
        return;
      end
      
      
      
      % easy way to get all data
      if isequal(varargin, {'ALL'})
        % if just 'ALL' is specified, the entire spectrogram is returned
        spec.(dType) = obj.mfObj.(dType);
        spec.Freqs = obj.Freqs;
        spec.Times = obj.Times;
        spec.WaveletParams = obj.WaveletParams;
        return;
      end
      
      % determine if the user wants to use samples or time points. Use time
      % points by default
      if any(strcmp(varargin,'SAMPLES'))
        timeBased = false;
      else
        timeBased = true;
      end
      
%       % select channels
%       if any(strcmp(varargin,'CHANINDICES'))
%         chanInds = varargin{find(strcmp(varargin,'CHANINDICES'))+1};
%         if any(chanInds<0) || any(chanInds>obj.NumChans)
%           error('Channel indices out of range');
%         end
%       else
%         chanInds = 1:obj.NumChans;
%       end
      
%       % select channels by name
%       if any(strcmp(varargin,'CHANNAMES'))
%         chanNames = varargin{find(strcmp(varargin,'CHANNAMES'))+1};
%         
%         try
%           chanInds = cellfun(@(x)find(strcmp(x,obj.ChannelNames)),chanNames);
%         catch
%           error('Channel names invalid');
%         end
%         
%         if any(chanInds<0) || any(chanInds>obj.NumChans)
%           error('Channel indices out of range');
%         end
%       else
%         chanInds = 1:obj.NumChans;
%       end
      
      % select frequencies
      if any(strcmp(varargin,'FREQINDICES'))
        % select particular frequencies by their index
        freqInds = varargin{find(strcmp(varargin,'FREQINDICES'))+1};
        if any(freqInds<0) || any(freqInds>obj.NumFreqs)
          error('Frequency indices out of range');
        end
      elseif any(strcmp(varargin,'FREQWINDOW'))
        % select a band of frequencies, 1x2 array with [freqLow freqHigh]
        freqBand = varargin{find(strcmp(varargin,'FREQWINDOW'))+1};
        freqInds = find((obj.Freqs >= freqBand(1))&(obj.Freqs <= freqBand(2)));
      else
        freqInds = 1:obj.NumFreqs;
      end
      
      % select time points
      if any(strcmp(varargin, {'SEGS'}))
        % select particular time periods, Nx2 array, first colunn is
        % segment starts, second column is segment ends.
        tWindows = varargin{find(strcmp(varargin,{'SEGS'}))+1};
        if ~isnumeric(tWindows) || (size(tWindows,2)~=2)
          error('Time windows are not properly formatted');
        end
        
        for j = 1:size(tWindows,1)
          if timeBased
            startInd = FindTimeIndex(obj,tWindows(j,1));
            endInd = FindTimeIndex(obj,tWindows(j,2));
          else
            startInd = tWindows(j,1);
            endInd = tWindows(j,2);
          end
          spec(j).(dType) = obj.mfObj.(dType)(freqInds,startInd:endInd);
          spec(j).Freqs = obj.Freqs(1,freqInds);
          spec(j).Times = obj.mfObj.sampleTimes(1,startInd:endInd);
          spec(j).WaveletParams = obj.WaveletParams;
        end
      elseif any(strcmp(varargin, {'PTS'}))
        % select particular time points, Nx1 array
        tPts = varargin{find(strcmp(varargin,{'PTS'}))+1};
        currInds = [];
        tList = obj.mfObj.sampleTimes;
        if timeBased
          if any(tPts > obj.mfObj.sampleTimes(1,end)) || any(tPts < obj.mfObj.sampleTimes(1,1))
            error('Time points outside of sampled times');
          end
          
          for j = 1:length(tPts)
            currInds(j) = FindTimeIndex(obj,tPts(j),tList);
          end
        else
          if any(tPts < 0) || any(tPts > obj.NumSamps)
            error('Time points outside of sampled times');
          end
          currInds = tPts;
        end
        spec.Freqs = obj.Freqs(1,freqInds);
        
        % preallocate spectrum
        spec.(dType) = zeros(length(spec.Freqs),length(tPts));
        spec.Times = zeros(1,length(tPts));
        spec.WaveletParams = obj.WaveletParams;
        for j = 1:length(currInds)
          spec.(dType)(:,j) = obj.mfObj.(dType)(freqInds,currInds(j));
          spec.Times(j) = tList(1,currInds(j));
        end
      elseif any(strcmp(varargin, {'WINDOWS'}))
        % select particular windows of data with a reference time and fixed
        % length of time (i.e. samples) before and after. First entry is
        % the length of time to sample before the reference point, the
        % second is the length of time to sample after the reference point,
        % and the remaining entries are the reference point times.
        tWinParams = varargin{find(strcmp(varargin,{'WINDOWS'}))+1};
        
        refTimes = tWinParams(3:end);
        if timeBased
          sampsBefore = floor(tWinParams(1)/obj.TimeStep);
          sampsAfter = ceil(tWinParams(2)/obj.TimeStep);
          if ((FindTimeIndex(obj,min(refTimes))-sampsBefore) <= 0)
            error('Window exceeds beginning of spectrum');
          elseif ((FindTimeIndex(obj,max(refTimes))+sampsAfter) > obj.NumSamps)
            error('Window exceeds end of spectrum');
          end
        else
          sampsBefore = tWinParams(1);
          sampsAfter = tWinParams(2);
          if (min(refTimes)-sampsBefore) <= 0
            error('Window exceeds beginning of spectrum');
          elseif (max(refTimes)-sampsAfter) > obj.NumSamps
            error('Window exceeds end of spectrum');
          end
        end
        
        for j = 1:length(refTimes)
          if timeBased
            currInd = FindTimeIndex(obj,refTimes(j));
          else
            currInd = refTimes(j);
          end
          spec(j).Freqs = obj.Freqs(1,freqInds);
          spec(j).Times = obj.mfObj.sampleTimes(1,(currInd-sampsBefore):(currInd+sampsAfter));
          spec(j).(dType) = obj.mfObj.(dType)(freqInds,(currInd-sampsBefore):(currInd+sampsAfter));
          spec(j).WaveletParams = obj.WaveletParams;
        end
      else
        spec.Freqs = obj.Freqs(1,freqInds);
        spec.Times = obj.mfObj.sampleTimes;
        spec.(dType) = obj.mfObj.(dType)(freqInds,:);
        spec.WaveletParams = obj.WaveletParams;
      end
    end
    
  end
end