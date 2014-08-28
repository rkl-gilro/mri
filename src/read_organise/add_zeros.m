function [num_zeros zeros_str] = add_zeros(num, digit_length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the number of zeros to be added given the number and the length
%% of the digit
%%
%% Inputs:  1. num          -> (scalar) number
%%          2. digit_length -> (scalar) the number of the digit
%% Outputs: 1. num_zeros    -> (scalar) number of added zeros
%%          2. zeros_str    -> (string) number of zeros
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expon = 1:digit_length;
array = 10 .^ expon;

for i = 1:length(array)
    if num < array(i)
        num_zeros = digit_length - i;
        
        tmp = [];
        for j = 1:num_zeros
            tmp = strcat(tmp, '0');
        end
        
        zeros_str = tmp;
        break;
    end
end