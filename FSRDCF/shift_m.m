function [ shift ] = shift_m( m )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    shift = eye(m+1);
    shift =flipud(shift);
    shift(1,1) = 1;
    shift = shift(1:m,1:m);
end

