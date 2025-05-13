function [gd] = xyz2gd(xyz)
%지구중심좌표계(ECEF)를 측지좌표계(Geodetic)로 변환

gd = NaN(size(xyz,1), 3);

%WGS-84 제원
a=6378137.0;
f = 1/298.257223563;
b = a*(1. - f);
aSQ = a^2;
bSQ = b^2;
eSQ = (aSQ - bSQ)/aSQ;

for i=1:size(xyz,1)
    x = xyz(i,1); y = xyz(i,2); z = xyz(i,3);

    %Computation of longitude
    longi = atan2(y,x)*180/pi;
    
    if longi > 180.
        longi = longi - 360.;
    elseif longi < -180.
        longi = longi + 360.;
    end
    
    %Iterative computaions of latitude and height
    p = sqrt(x^2 + y^2);
    q = 0;
    
    Phi0 = atan2(z*inv(1-eSQ), p);
    
    while (q ~= 1)
        NO = aSQ/sqrt(aSQ*(cos(Phi0))^2 + bSQ*(sin(Phi0))^2);
        h = p/cos(Phi0)-NO;
        Phi = atan2(z*inv(1 - eSQ*(NO/(NO + h))), p);
    
        if abs(Phi - Phi0) <= 1e-13
            break
        else
            Phi0 = Phi;
        end
    end
    
    lat = Phi*180/pi;
    
    gd(i,:) = [lat, longi, h];
end