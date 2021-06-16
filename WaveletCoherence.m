function wcoh = WaveletCoherence(spec1, spec2, varargin)
% spec1 - a spectrogram structure or WaveletObj, use WaveletSpec to generate spectrogram
% spec2 - a spectrogram structure or WaveletObj, use WaveletSpec to generate spectrogram

% Optional arguments
% SAVE - followed by a string, saves the file to the filename given by
% the string

if isa(spec1,'WaveletObj')
    isObj1 = true;
    spec1Size = {[spec1.NumFreqs spec1.NumSamps]};
else
    isObj1 = false;
    spec1Size = arrayfun(@(x)size(x.Spectrum),spec1,'UniformOutput',false);
end

if isa(spec2,'WaveletObj')
    isObj2 = true;
    spec2Size = {[spec2.NumFreqs spec2.NumSamps]};
else
    isObj2 = false;
    spec2Size = arrayfun(@(x)size(x.Spectrum),spec2,'UniformOutput',false);
end

% check for consistency between spectrograms
if ~isequal(spec1Size, spec2Size)
    error('Spectrograms are of unequal sizes');
elseif ~isequal({spec1.Freqs},{spec2.Freqs})
    error('Spectrograms do not contain matching frequencies');
elseif ~isequal({spec1.WaveletParams}, {spec2.WaveletParams})
    error('Spectrograms have different wavelet parameters');
end

if any(strcmp(varargin,'SAVE'))
    if (numel(spec1Size) > 1) || (numel(spec2Size) > 1)
        error('Cannot save coherence spectrograms for arrays');
    end
    fileName = varargin{find(strcmp(varargin,'SAVE'))+1};
    if ~strcmp(fileName((end-3):end),'.wls')
        fileName = [fileName '.wls'];
    end
    saveYes = true;
else
    saveYes = false;
end



if saveYes
    if length(spec1Size) > 1
        error('Can only save when a single time window is used');
    end
    Coherence = zeros(1,spec1Size{1}(2));
    Phase = zeros(1,spec1Size{1}(2));
    waveletParams = spec1.WaveletParams;
    freqs = spec1.Freqs;
    sampleTimes = spec1.Times;
    save(fileName,'sampleTimes', 'freqs', 'Coherence', 'Phase', ...
        'waveletParams', '-v7.3');
    clear sampleTimes coherence phase waveletParams freqs;
    mfObj = matfile(fileName, 'Writable',true);
    wcoh = fileName;
end

for k = 1:numel(spec1Size)
    
    % determine wavelet type for generating smoothing window
    wType = spec1(k).WaveletParams{1};
    wParams = spec1(k).WaveletParams{2};
    
    for j = 1:length(spec1(k).Freqs)
        currFreq = spec1(k).Freqs(j);
        if isObj1
            currSpec1 = spec1.LoadSpectrum('FREQINDICES',j);
            currSpec1 = currSpec1.Spectrum;
        else
            currSpec1 = spec1(k).Spectrum(j,:);
        end
        
        if isObj2
            currSpec2 = spec2.LoadSpectrum('FREQINDICES',j);
            currSpec2 = currSpec2.Spectrum;
        else
            currSpec2 = spec2(k).Spectrum(j,:);
        end
        
        
        tStep = median(diff(spec1(k).Times));
        switch wType
            case 'morlet'
                sigma = wParams(1)/(2*pi*currFreq);
                amp = 1/(sigma*sqrt(2*pi));
                tSupp = (-wParams(1)/currFreq):tStep:(wParams(1)/currFreq);
                cmwl = amp*exp(-(tSupp.^2)./(2*(sigma^2))).*exp(2*1i*pi*currFreq*tSupp);
                w1 = real(cmwl);
                w2 = imag(cmwl);
            otherwise
                error('No valid wavelet');
        end
        
        % Smooth spectrograms
        
        sigma = 7/(2*pi*currFreq);
        amp = 1/(sigma*sqrt(2*pi));
        tSupp = (-7/currFreq):tStep:(7/currFreq);
        smKern = amp*exp(-(tSupp.^2)./(2*(sigma^2)));
        smKern = smKern/sum(smKern);
        s1SmMag = conv(abs(currSpec1).^2,smKern,'same');
        s2SmMag = conv(abs(currSpec2).^2,smKern,'same');
        
        % cross of spectrums
        xS1S2 = currSpec1 .* conj(currSpec2); %spec1.s.*conj(spec2.s);
        xS1S2Sm = conv(real(xS1S2),smKern,'same') + ...
            1i*conv(imag(xS1S2),smKern,'same');
        
        % calculate phase and coherence
        xSpec = xS1S2Sm./(sqrt(s1SmMag).*sqrt(s2SmMag));
        currCoh = (abs(xS1S2Sm).^2)./(s1SmMag.*s2SmMag);
        currPhase = angle(xSpec);
        
        % store results
        if saveYes
            mfObj.Coherence(j,:) = currCoh;
            mfObj.Phase(j,:) = currPhase;
        else
            wcoh(k).Times = spec1(k).Times;
            wcoh(k).Freqs = spec1(k).Freqs;
            wcoh(k).WaveletParams = spec1(k).WaveletParams;
            wcoh(k).Coherence(j,:) = currCoh;
            wcoh(k).Phase(j,:) = currPhase;
        end
    end
end
end