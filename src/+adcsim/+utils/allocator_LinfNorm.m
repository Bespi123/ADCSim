function T_winf = allocator_LinfNorm(T_w_L2norm)
%Function that calculates an allocator based on the L2-norm
% Calculate the pseudoinverse matrix for the W with the 
% Moore-Penrose Pseudoinvers
    alpha = -1*(min(T_w_L2norm)+max(T_w_L2norm))/2;
    T_winf = T_w_L2norm + alpha*[1,1,1,1]';
end
