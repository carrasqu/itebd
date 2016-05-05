function [ energy ] = tedben( LB,lambB,LA,lambA,hl,chimax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lb2=lambB*lambB;
energy=scon({lb2,LA,lambA,LB,lb2,LB,lambA,LA,hl},{[1 8] ,[1 11 2] ,[2 3] ,[3 12 4] ,[ 4 5], [6 10 5], [ 7 6 ] ,[ 8 9 7], [9 10 11 12]});
norm=scon({lb2,LA,lambA,LB,lb2,LB,lambA,LA},{[1 8] ,[1 9 2] ,[2 3] ,[3 10 4] ,[ 4 5], [6 10 5], [ 7 6 ] ,[ 8 9 7]});

end

