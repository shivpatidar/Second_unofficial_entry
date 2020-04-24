function V = centtempcorren(x,L,sigma)
% CENTTEMPCORREN calculates the centered temporal correntropy of a time series
% Input: - x: time series
%		 - L: number of lags to consider
%		 - sigma: Gaussian kernel size
% Output: V, centered temporal correntropy matrix of size LxL
%
% Comments: The code uses Incomplete Cholesky Decomposition by Sohan Seth.
%
% Author: Steven Van Vaerenbergh (steven@gtas.dicom.unican.es) Date: 17.11.2008
%
% USAGE: V = centtempcorren(x,L,sigma)

x = x(:);
N = length(x);	% number of data points

% calculate the cholesky decomposition
G =  incompleteCholesky(x,sigma);

% fill the correntropy matrix
V = zeros(L);
num = N-L+1;
num2 = num^2;
% inum2 = 1/num2;
twoSigmaSquare = 2*sigma^2;
% snum1 = sqrt(num-1);
z = zeros(N,1);
for i=1:L,			% diagonal indicator in lag index matrix
% 	if (mod(i,floor(L/10))==1), fprintf(1,'.'); end

	for j=0:L-i,	% index on diagonal in lag index matrix
		ind1 = 1+j:1+j+num-1;
		ind2 = i+j:i+j+num-1;
		if (j==0)
			% initialize first submatrix for the "lag diagonal"
			v = num*sum(exp(-(x(ind1)-x(ind2)).^2/twoSigmaSquare));
			o1_ini = z; o1_ini(ind1) = 1;
			o2_ini = z; o2_ini(ind2) = 1;
			v_mean = (o1_ini'*G)*(G'*o2_ini);
			V(ind1(1),ind2(1)) = v-v_mean;

			o1 = z;
			o2 = z;
			o1(ind1(2:end)) = 1;
			o2(ind2(2:end)) = 1;
		else
			% efficiently compute correntropy for next submatrix on "lag diagonal"
			Go1_nr = G(ind1(end),:)';
			Go2_oc = G(ind2(1)-1,:)';
			Go2_nc = G(ind2(end),:)';
			Go1_or = G(ind1(1)-1,:)';
			V(ind1(1),ind2(1)) = V(ind1(1)-1,ind2(1)-1)...
				+ (Go1_nr'*Go2_nc - Go1_or'*Go2_oc)*(num-1)...
				+ o2'*G*(Go1_or-Go1_nr)...
				+ o1'*G*(Go2_oc-Go2_nc);
			
			% update indices for next iteration
			o1(ind1(1)) = 0;
			o1(ind1(end)) = 1;
			o2(ind2(1)) = 0;
			o2(ind2(end)) = 1;
			
		end
		
	end
end
% fprintf(1,'\n');

V = V + V';
V = V - diag(diag(V))/2;
V = V/sqrt(2*pi)/sigma/num2;

% ~ Done