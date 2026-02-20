function rotX = quatRotation(q, x)
    % QUATROTATION Rotates a 3D vector using an attitude quaternion.
    %
    %   This function applies Hamilton's rotation operator (sandwich 
    %   product) to transform a vector from one reference frame to another.
    %   Mathematical formula: x' = q * [0, x] * q_conjugate
    %
    %   Inputs:  
    %       q : Rotation quaternion [qw, qx, qy, qz] (1x4 or 4x1).
    %       x : 3D vector to rotate [x, y, z] (1x3 or 3x1).
    %
    %   Output:
    %       rotX : Rotated 3D vector (1x3).
    
    %%% Import adsim utilities
    import adcsim.utils.*

    % 1. Convert the 3D vector into a "pure quaternion"
    % A 0 is added to the scalar part to enable quaternion operations.
    qx = [0, x(1), x(2), x(3)];
    
    % 2. Ensure the quaternion 'q' is a row vector (1x4)
    [a, b] = size(q);
    if a == 4 && b == 1
        q = q'; % Transpose from column to row
    end

    % 3. Apply the full rotation (Sandwich product)
    % First multiplies q * qx, then multiplies the result by the conjugate of q.
    qrotX = quaternProd(quaternProd(q, qx), quaternConj(q));
    
    % 4. Extract the resulting vector
    % The scalar part (index 1, which will be very close to 0) is discarded, 
    % keeping only the spatial components x, y, and z.
    rotX = qrotX(2:4);
end