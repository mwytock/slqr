function [A,B,W,Q,R,K_mask] = problem_wac(filename)
load(['~/energy/control/wac_pst/' filename '.mat']);
W = B*B';
