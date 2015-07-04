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
        % x should be a real number (not array) within the domain of the spline
        function r = eval(obj, x)
            r = Splinter.getInstance().call('eval', obj.Handle, x, length(x));
            %r = calllib(Splinter.getInstance.get_alias(), obj.Eval_function, obj.Handle, x, length(x));
        end
        
        % Evaluate the Jacobian at x
        % x should be a real number (not array) within the domain of the spline
        function r = evalJacobian(obj, x)
            temp = Splinter.getInstance().call('eval_jacobian', obj.Handle, x, length(x));
            reshape(temp, 1, length(x))
            r = temp.value;
        end
        
        % Evaluate the Hessian at x
        % x should be a real number (not array) within the domain of the spline
        function r = evalHessian(obj, x)
            temp = Splinter.getInstance().call('eval_hessian', obj.Handle, x, length(x));
            reshape(temp, length(x), length(x))
            r = temp.value;
        end

        function r = getNumVariables(obj)
            r = Splinter.getInstance().call('get_num_variables', obj.Handle);
        end
        
        function save(obj, fileName)
            Splinter.getInstance().call('save', obj.Handle, fileName);
        end
        
        function load(obj, fileName)
            Splinter.getInstance().call('load', obj.Handle, fileName);
        end

        % Destructor. Deletes the internal Approximant object.
        function delete(obj)
            if(obj.Handle ~= -1)
                Splinter.getInstance().call('delete_approximant', obj.Handle);
            end
        end
    end
end

