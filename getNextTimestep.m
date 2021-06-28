function nextTimeStep = getNextTimestep(sita, Lp, Lo, x, kbT) 
       constantCoeff = kbT/(Lp*Lo);
       change = 0.5/((1 - x/Lo)^3) + 1;
        gradient = abs(constantCoeff*change);
        nextTimeStep = (2*0.01 * sita) / gradient;
        
       

        