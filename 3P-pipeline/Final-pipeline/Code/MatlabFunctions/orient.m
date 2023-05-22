function [objOrientation, angleObj2Major] = orient(majOrientation, objAngle)
     % This function takes the orientation of the major axis of the 
     % cell (majOrientation, range [-90 90] degrees) and the angle 
     % formed by a line drawn from the centroid of the cell to the 
     % centroid of an object (objAngle, [0 360] degrees), and 
     % converts the objAngle to objOrientation ([-90 90] degrees,
     % and finds the value of the acute angle formed by the
     % intersection of the line with the major axis.

    if isnan(objAngle)
        objOrientation = NaN;
        angleObj2Major = NaN;
    else
        if objAngle > 270
            objOrientation = objAngle - 360;
        elseif objAngle > 90 
            objOrientation = objAngle - 180;
        else
            objOrientation = objAngle;
        end

        diffAngle = abs(objOrientation - majOrientation);
        if diffAngle <= 90
            angleObj2Major = diffAngle;
        elseif diffAngle > 180 
            error('diffAngle > 180')
        elseif diffAngle > 90
            angleObj2Major = 180 - diffAngle;       
        end  
        if angleObj2Major < 0
            error('angleObj2Major is negative')
        end
    end
end
