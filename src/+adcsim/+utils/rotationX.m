function Rx = rotationX(theta)
    %ROTATIONX Generates a 3x3 rotation matrix around the X-axis.
    %
    %   Rx = ROTATIONX(THETA) returns a 3x3 rotation matrix that represents
    %   a counter-clockwise rotation by an angle THETA (in radians) around
    %   the X-axis.
    %
    %   Inputs:
    %     theta - The angle of rotation in radians.
    %
    %   Output:
    %     Rx    - A 3x3 rotation matrix.
    %
    
    % Construct the rotation matrix
    Rx = [1,           0,           0;
          0,  cos(theta), -sin(theta);
          0,  sin(theta),  cos(theta)];
end