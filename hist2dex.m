
mYX = rand(100,2); 
vXEdge = linspace(0,1,10); 
vYEdge = linspace(0,1,20); 
mHist2d = hist2d(mYX,vYEdge,vXEdge); 
  
nXBins = length(vXEdge); 
nYBins = length(vYEdge); 
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins)); 
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins)); 
pcolor(vXLabel, vYLabel,mHist2d); colorbar