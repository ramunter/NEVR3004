function[yy]= interpolateBetweenPoints(xx,ind0,ind1,y0,y1,angleFreq)               %ind0 and ind1 is the two closest indices of angleData 

 

% if this looks stange, google "linear interpolation" in wikipedia

yd=y1-y0;
if(isnan(y0) || isnan(y1))   %if firingTime is close to a nan angle, we do not interpolate 
    yy=nan;
else 
    if(yd>pi)               % 
        y1=y1-2*pi;
    elseif(yd<-pi)
        y1=y1+2*pi;
    end
   
    x0=1000*ind0/angleFreq;            
    x1=1000*ind1/angleFreq;
    
    yy=y0 + (xx-x0)*(y1-y0)/(x1-x0); %linear interpolation
end
   
yy=mod(yy,2*pi);
 
end


