% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef Utilities
    methods(Static)
        % Returns a (column-major) (double) C array with the same contents as mat
        % Note: The returned c array is actually using the same underlying
        % data as MatLab, so be careful not to modify it in C (unless you
        % know what you're doing, in which case you wouldn't have bothered
        % about this warning anyway).
        function r = matrixToCArray(mat)
            libP = libpointer('doublePtr', mat);
            reshape(libP, 1, numel(mat));
            r = libP;
        end
        
        % Returns a MatLab matrix of size numRowsxnumCols with the data
        % from the C array of doubles cArray.
        % The underlying data is probably the same as is the case with
        % matrixToCArray, but we have not researched it, so don't take our
        % word for it.
        function r = cArrayToMatrix(cArray, numRows, numCols)
            reshape(cArray, numRows, numCols);
            r = cArray.value;
        end
        
        % Returns a MatLab matrix of size numRowsxnumColsxnumPages with the
        % data from the C array of doubles cArray.
        % The underlying data is probably the same as is the case with
        % matrixToCArray, but we have not researched it, so don't take our
        % word for it.
        function r = cArrayTo3dMatrix(cArray, numRows, numCols, numPages)
            reshape(cArray, numRows, numCols*numPages);
            r = reshape(cArray.value, numRows, numCols, numPages);
        end
    end
end