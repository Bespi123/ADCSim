function ts = calculateSettlementTime(e, t, tol)
% calculateSettlementTime computes the settlement time for a given signal `e`.
% The function checks when the signal `e` first enters and stays within the defined 
% tolerance range for the entire signal and returns the corresponding time at that 
% point. If the signal never stays within the range, it returns NaN.
%
% Inputs:
%   e   - A vector representing the signal values (e.g., angular error or any monitored variable)
%   t   - A vector representing the time instances corresponding to each value in `e`.
%   tol - A scalar tolerance value defining the acceptable range for the signal `e`.
%
% Output:
%   ts  - The settlement time, i.e., the time at which the signal `e` stays within 
%         the tolerance range [-tol, tol]. If the signal never settles, `ts` will be NaN.
    
    %%% Set bounds based on the given tolerance (lower and upper bounds)
    bounds = [-tol, tol];  % The bounds are set as [-tol, tol]
    
    %%% Find the indices where the signal `e` is outside the bounds
    % The expression `~(e < bounds(2) & e > bounds(1))` finds values outside the range [-tol, tol]
    [x,~] = find(~(e < bounds(2) & e > bounds(1))); 
    
    %%% Check if the last value of `x` corresponds to the final time index
    if (x(end) == length(t))  % If the last index where the signal is out of range is the final time
        ts = NaN;  % If the signal never settles within the tolerance, return NaN
    else
        %%% If the signal settles within the range, return the corresponding time
        ts = t(max(x));  % Return the time at the last instance where the signal was out of range
    end
end
