function rotX = quatRotation(q, x)
    % QUATROTATION Vectorized rotation of 3D vectors using quaternions.
    %
    %   This function rotates a set of 3D vectors by a set of quaternions
    %   using the efficient vector-form of the sandwich product.
    %
    %   INPUTS:
    %       q - (Nx4) Matrix of quaternions [qw, qx, qy, qz]
    %       x - (Nx3) Matrix of vectors to rotate [vx, vy, vz]
    %
    %   OUTPUT:
    %       rotX - (Nx3) Matrix of rotated vectors.
    %
    %   METHOD:
    %       Uses the formula: v' = v + 2*s*(r x v) + 2*(r x (r x v))
    %       where q = [s, r] (s is scalar, r is vector part).

    % 1. Ensure inputs are matrices of size (Nx4) and (Nx3)
    if size(q, 2) ~= 4, q = q'; end
    if size(x, 2) ~= 3, x = x'; end
    
    % 2. Extract scalar (s) and vector (r) parts of the quaternion
    s = q(:, 1);             % (Nx1)
    r = q(:, 2:4);           % (Nx3)
    
    % 3. Calculate intermediate cross products
    % cross(r, x, 2) performs the cross product row-wise
    r_cross_x = cross(r, x, 2);                      % (Nx3)
    r_cross_r_cross_x = cross(r, r_cross_x, 2);      % (Nx3)
    
    % 4. Apply the rotation formula
    % rotX = x + 2*s*(r x x) + 2*(r x (r x x))
    rotX = x + 2 .* s .* r_cross_x + 2 .* r_cross_r_cross_x;
end