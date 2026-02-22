function Rz = rotationZ(theta)
    %ROTATIONZ Generates a 3x3 rotation matrix around the Z-axis.
    %
    %   Rz = ROTATIONZ(THETA) returns a 3x3 rotation matrix that represents
    %   a counter-clockwise rotation by an angle THETA (in radians) around
    %   the Z-axis.
    %
    %   Inputs:
    %     theta - The angle of rotation in radians.
    %
    %   Output:
    %     Rz    - A 3x3 rotation matrix.
    %
    
    % Construct the rotation matrix
    Rz = [cos(theta), -sin(theta), 0;
          sin(theta),  cos(theta), 0;
          0,           0,          1];
    
end