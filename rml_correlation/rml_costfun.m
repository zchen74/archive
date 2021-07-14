%---------------------------------------------------------------------------%
% RML cost function 							    %
%---------------------------------------------------------------------------%

% theta(1) = the variance
% theta(2) = temporal covariance
% theta(3) = spatial covariance


function [ L ] = rml_costfun(theta, z,X,D,E,ntimes)

	% disp('Current parameter guess');
	disp(num2str(theta));
	% tic;

	% Reset the decorrelation time to 20 days if it exceeds 20 days.
	% A decorrelation time >20 days will cause D to become non-invertible.
	if theta(2) > 20; theta(2) = 20; end;


%-------------------------------------------------------------%
% Multiply the distance matrices by the covariance functions  %
%-------------------------------------------------------------%

	% Temporal covariance: Spherical covariance model
%	D1 = 1 - 1.5 .* (D   ./theta(2))  + 0.5 .* (D.^3   ./ theta(2).^3);
%	D1(D   > theta(2)) = 0;

	% Spatial covariance: spherical model
%	E1 = 1 - 1.5 .* (E   ./theta(3))  + 0.5 .* (E.^3   ./ theta(3).^3);
%	E1(E   > theta(3)) = 0;

        D1 = exp(-1 .* D ./theta(2));
        D1(D > 3.*theta(2)) = 0;
        E1 = exp(-1 .* E ./theta(3));
        E1(E > 3.*theta(3)) = 0;  

%---------------------------------%
% Calculate the log determinant   %
%---------------------------------%

	% L1 = 0.5.*log(det(csigma));

	% More computationally efficient formulation of L1:
	% Really slow!
	%L1 = 0.5 .* 2.*sum(log(diag(chol(theta(1) .* kron(D1,E1)))));
	
	% Much faster option (but could still be faster)
	L1 = 0.5 .* 2.*sum(log(kron(diag(chol(D1)),diag(chol(theta(1).*E1)))));
	
	% Probably even faster
        % L1 = 0.5 .* 2.*sum(log(diag(kron(diag(chol(D1)),diag(chol(theta(1).*E1))))));


%------------------------------------------------%
% Calculate second portion of the cost funciton  %
%------------------------------------------------%

	% Here, I need to calculate X'QX

	% Calculate QX (We'll call it B)
	% B = Q * X
	B = [];
	
	m1 = size(X,1) ./ ntimes; 
	for j = 1:size(D1,1);
	B1 = zeros(m1,size(X,2));
		for i = 1:size(D1,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		B1 = B1 + X(sel,:) .* D1(j,i) .* theta(1).*theta(1);
		end; % End of i loop
	temp = E1 * B1;
	B = [B; temp];
	end; % End of i loop
	clear B1 temp;

	% Calculate the 2nd portion of the cost function
	L2 = 0.5.*log(det(X.'*B));


%-----------------------------------------------------%
% Calculate the third component of the cost function  %
%-----------------------------------------------------%

	% Calculate G * z
	% G * z = A - B * inv( X' * B) * B' * z

	% Take the inverse of D and E	
	Dinv = inv(D1);
	Einv = inv(E1);


	%---------------%
	% Calculate A   %
	%---------------%
	
	% A = inv(Q) * z
	A = [];
	
	for j = 1:size(Dinv,1);
	A1 = zeros(m1,1);
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		A1 = A1 + z(sel) .* Dinv(j,i) .* (1./theta(1).*(1./theta(1)));	
		end; % End of i loop
	temp = Einv * A1;
	A = [A; temp];
	end; % End of i loop
	clear A1 temp;
	
	
	%---------------%
	% Calculate B   %
	%---------------%
	
	% B = inv(Q) * X
	B = [];
	
	for j = 1:size(Dinv,1);
	B1 = zeros(m1,size(X,2));
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		B1 = B1 + X(sel,:) .* Dinv(j,i) .* (1./theta(1).*(1./theta(1)));
		end; % End of i loop
	temp = Einv * B1;
	B = [B; temp];
	end; % End of i loop
	clear B1 temp;
	
	
	%------------------------%
	% Calculate z' * G * z   %
	%------------------------%
	
	Gz = A - B * ((X' * B) \ (B' * z));
	
	L3 = 0.5 .* z' * Gz;
	clear Gz;


%--------------------------------------%
% Assemble the final cost function     %
%--------------------------------------%

	L = L1 + L2 + L3;

	% disp('Cost function value');
	% disp(num2str(L));

	% disp('Time elapsed');
	% disp(toc);
