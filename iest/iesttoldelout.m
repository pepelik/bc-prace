function [x,y] = iesttoldelout(x, y, deltas, n)
%BEGINDOC=================================================================
% .Author.
%
%  Petra Pelikanova
%
%-------------------------------------------------------------------------
% .Description.
%
%   Function deletes from dataset n outliers (the measurement which corespond
%   the biggest delta parameters).
%
%-------------------------------------------------------------------------
% .Input parameters.
%
%   x ... input data vector
%   y ... vector of output data
%   deltas ... delta parameters from iesttol()
%   n ... number of outliers
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   [x, y] ... data without n outliers
%
%------------------------------------------------------------------------
% .Implementation details.
%
%   See Hladik, Milan, and Michal Cerny. "Interval regression by tolerance
%   analysis approach." Fuzzy Sets and Systems 193 (2012): 85-107.
%
%
%------------------------------------------------------------------------
%  .Licence.
%
%   Copyright (C) 2015  Charles University in Prague, Czech Republic
%
%   LIME 1.0 is free for private use and for purely academic purposes.
%   It would be very kind from the future user of LIME 1.0 to give
%   reference that this software package has been developed
%   at Charles University, Czech Republic.
%
%   For any other use of LIME 1.0 a license is required.
%
%   THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, WITHOUT LIMITATIONS, THE IMPLIED WARRANTIES
%   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%
%-------------------------------------------------------------------------
% .History.
%
%   2017-08-08    first version
%
%------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC===================================================================

y_len = length(y);
sorted_deltas = sort(deltas);
max_delta = sorted_deltas(y_len-n+1);
for j = flip(1:y_len)
	if(deltas(j)>=max_delta && n > 0)
		y = [y(1:j-1); y(j+1:end)];
		x = [x(1:j-1,:); x(j+1:end,:)];
	end
end
