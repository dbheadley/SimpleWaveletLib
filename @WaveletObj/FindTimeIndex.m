%% uses a binary search to find the closest value
function tInd = FindTimeIndex(obj, tVal, varargin)
 
  if ~isempty(varargin)
    tList = varargin{1};
  else
    tList = obj.mfObj.sampleTimes;
  end
  
  numTVals = obj.NumSamps;
  window = [1 numTVals];
  
  while (window(2)-window(1))>1
    compInd = window(1) + round((window(2)-window(1))/2);
    compVal = tList(1,compInd);
    
    if tVal == compVal
      tInd = compInd;
      return;
    elseif tVal < compVal
      window(2) = compInd-1;
    else
      window(1) = compInd+1;
    end
  end

  
  if abs(tList(1,window(1))-tVal) > ...
      abs(tList(1,window(2))-tVal)
    tInd = window(2);
  else
    tInd = window(1);
  end
  
  if tInd ~= 1
    if abs(tList(1,tInd-1)-tVal) < ...
        abs(tList(1,tInd)-tVal)
      tInd = tInd - 1;
    end
  end
  
  if abs(tList(1,tInd)-tVal) > ...
      abs(tList(1,1 + round((numTVals-1)/2)) - tVal)
      tInd = 1 + round((numTVals-1)/2);
  end
end