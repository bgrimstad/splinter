% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef (Abstract = true) Approximant < handle
    properties (Abstract = true, Access = protected)
        Handle
    end

    methods
        % Evaluate the spline at x
        % Supports batch evaluation:
        % x = [0 0; 1 1] will evaluate the approximant in both [0 0] and
        % [1 1], and return a nx1 matrix with the values, where n is the
        % number of rows in x (or points you evaluated it in).
        function r = eval(obj, x)
            libP = Utilities.matrixToCArray(x);
            numPoints = numel(x) / obj.getNumVariables();
            temp = Splinter.getInstance().call('eval_col_major', obj.Handle, libP, numel(libP.value));
            r = Utilities.cArrayToMatrix(temp, numPoints, 1);
        end
        
        % Evaluate the Jacobian at x
        function r = evalJacobian(obj, x)
            libP = Utilities.matrixToCArray(x);
            numPoints = numel(x) / obj.getNumVariables();
            temp = Splinter.getInstance().call('eval_jacobian_col_major', obj.Handle, libP, numel(libP.value));
            r = Utilities.cArrayToMatrix(temp, numPoints, obj.getNumVariables());
        end
        
        % Evaluate the Hessian at x
        function r = evalHessian(obj, x)
            libP = Utilities.matrixToCArray(x);
            numPoints = numel(x) / obj.getNumVariables();
            temp = Splinter.getInstance().call('eval_hessian_col_major', obj.Handle, libP, numel(libP.value));
            r = Utilities.cArrayTo3dMatrix(temp, obj.getNumVariables(), obj.getNumVariables(), numPoints);
        end

        function r = getNumVariables(obj)
            r = Splinter.getInstance().call('get_num_variables', obj.Handle);
        end
        
        function save(obj, fileName)
            Splinter.getInstance().call('save', obj.Handle, fileName);
        end

        % Destructor. Deletes the internal Approximant object.
        function delete(obj)
            if(obj.Handle ~= -1)
                Splinter.getInstance().call('delete_approximant', obj.Handle);
            end
        end
    end
end