function X = lyap_fast(U, D, invU, Q);
X = real(-U*(invU*Q*invU.' ./ (D*ones(size(D)) + ones(size(D))*D))*U.');
