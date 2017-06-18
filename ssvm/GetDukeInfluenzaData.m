function [ X labels Y labelsH1N1] = GetDukeInfluenzaData( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clc; clear all; close all; randn('state',0); rand('state',0)
load Fluz.mat; load H1N1.mat;

ndx1 = find(Fluz.time==14 | Fluz.time==15 | Fluz.time==16);
ndx2 = find((H1N1.time==14 | H1N1.time==15 | H1N1.time==16) ...
    & H1N1.id~=1 & H1N1.id~=3 & H1N1.id~=5 & H1N1.id~=10 & H1N1.id~=18);
Train.data = Fluz.data1(ndx1,:)'; Train.label = Fluz.label(ndx1);
Test.data = H1N1.data1(ndx2,:)'; Test.label = H1N1.label(ndx2);
[p,n] = size(Train.data);
Train.mean = mean(Train.data,2);
X0 = Train.data - repmat(Train.mean,1,n);
Xscale = sqrt(sum(X0.*X0,2)); % Xscale = ones(p,1); %
X = [X0'./repmat(Xscale',n,1)]; o = Train.label;
labels  =o;


[p1,n1] = size(Test.data);
Test.mean = mean(Test.data,2);
Y0 = Test.data - repmat(Test.mean,1,n1);
Yscale = sqrt(sum(Y0.*Y0,2)); % Yscale = ones(p,1); % Yscale = Xscale; %
Y = [Y0'./repmat(Yscale',n1,1)];
labelsH1N1 = Test.label;

end

