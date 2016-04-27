function R = rotmat3(ang,comp)
%function R = rotmat3(ang,comp)
% where ang is the angle of the rotation and comp is either 1, 2, or 3 for
% the x, y, or z axis, indicating the axis of the rotation
% angles are in degrees and positive indicates CCW rotation
% if ang is a 3x1 vector and comp is 4, this will create a combined
% rotation matrix using ang(1) for the rotation about x, ang(2) for the
% rotation about y, and ang(3) for the rotation about z (really x prime).

if comp==1 && length(ang)==1
    R=[1   0          0;
       0   cosd(ang)  sind(ang);
       0  -sind(ang)  cosd(ang)];

elseif comp==2 && length(ang)==1
    R=[cosd(ang) 0 -sind(ang);
        0        1  0;
       sind(ang) 0  cosd(ang)];
    
elseif comp==3 && length(ang)==1
    R=[cosd(ang)  sind(ang)  0;
      -sind(ang)  cosd(ang)  0;
       0          0          1];
 % this is from Goldstein, Classical Mechanics, 3rd edition, p154  
elseif comp==4 && length(ang)==3
    cu=cosd(ang(1));
    ct=cosd(ang(2));
    cp=cosd(ang(3));
    su=sind(ang(1));
    st=sind(ang(2));
    sp=sind(ang(3));
    R=[cu*cp-ct*sp*su  -su*cp-ct*sp*cu  st*sp;
       cu*sp+ct*cp*su  -su*sp+ct*cp*cu -st*cp;
             st*su            st*cu     ct];

else
    error('comp must be 1, 2, 3, or 4')
end
