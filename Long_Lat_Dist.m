function [d_km]=Long_Lat_Dist(lonlat1,lonlat2)
% d_km: distance in km based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
radius=6371; % Earth radius in km
lon1=lonlat1(1)*pi/180;
lon2=lonlat2(1)*pi/180;
lat1=lonlat1(2)*pi/180;
lat2=lonlat2(2)*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c=2*atan2(sqrt(a),sqrt(1-a));
d_km=radius*c;    %Haversine distance
end