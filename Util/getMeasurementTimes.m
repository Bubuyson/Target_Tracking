function [MeasurementInterval] = getMeasurementTimes(tStart, tEnd, mean, var)

    thatTime = tStart;
    j = 0; eps = 1e-6;
    lastidtlim = max(0.1, mean - 1*sqrt(var)); flag = 0;
    while(1)
        if(tEnd - thatTime <= mean*1.1)
            if(tEnd - thatTime <= lastidtlim - eps), flag = 1; end
            break
        end
        R = mvnrnd_NT(mean, var);
        if((R < 0) || (thatTime + R > tEnd) ||(R > mean + 2*sqrt(var)))
            continue
        end
        j = j+1;
        Meas(j) = thatTime + R;     %#ok
        thatTime = thatTime + R;
    end
    if exist('Meas', 'var')
        MeasurementInterval = [tStart round(Meas, 4) tEnd];
        if flag, MeasurementInterval(end-1) = []; end
    else
        MeasurementInterval = [tStart tEnd];
    end

end