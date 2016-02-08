function [norm_ems, rev_n_ems] = amp_comp(m_values)
 %function [norm_ems, rev_n_ems] = amp_comp(m_values)
   %calculates normalized amplitude values, with both the maximum value as 1
   %and the minimum value as 1
     %Inputs:
       %m_values - 1 x numel(stations) vector of amplitude values at
         %selected time for phase
     %Outputs:
       %norm_ems - 1 x numel(stations) vector of the maximum value as 1 ->
         %normalized to the maximum value
       %rev_n_ems - 1 x numel(stations) vector of the minimum value as 1 ->
         %normalized to the minimum value (largest = smallest amplitudes)
ems = abs(m_values);
max_val = max(ems);
norm_ems = ems/max_val; %max value is 1
rev_norm_ems = 1-norm_ems;
max_val_norm = max(rev_norm_ems); 
rev_n_ems = rev_norm_ems/max_val_norm; %min value is 0 -> was max value