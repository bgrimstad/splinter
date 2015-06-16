classdef (Abstract = true) Approximant < handle
    properties (Abstract = true, Access = protected)
        Handle
        
        % Name of the functions in the back end, used with calllib
        Constructor_function
        Destructor_function
        Eval_function
        EvalJacobian_function
        EvalHessian_function
        Save_function
        Load_function
    end

    methods
        % Evaluate the spline at x
        % x should be a real number (not array) within the domain of the spline
        function r = eval(obj, x)
            r = Splinter.getInstance().call(obj.Eval_function, obj.Handle, x, length(x));
            %r = calllib(Splinter.getInstance.get_alias(), obj.Eval_function, obj.Handle, x, length(x));
        end
        
        % Evaluate the Jacobian at x
        % x should be a real number (not array) within the domain of the spline
        function r = evalJacobian(obj, x)
            temp = calllib(Splinter.getInstance().get_alias(), obj.EvalJacobian_function, obj.Handle, x, length(x));
            reshape(temp, 1, length(x))
            r = temp.value;
        end
        
        % Evaluate the Hessian at x
        % x should be a real number (not array) within the domain of the spline
        function r = evalHessian(obj, x)
            temp = calllib(Splinter.getInstance().get_alias(), obj.EvalHessian_function, obj.Handle, x, length(x));
            reshape(temp, length(x), length(x))
            r = temp.value;
        end
        
        function save(obj, fileName)
            calllib(Splinter.getInstance().get_alias(), obj.Save_function, obj.Handle, fileName);
        end

        % Destructor. Deletes the internal BSpline object.
        % Note that calling any of the methods after this,
        % WILL make MatLab crash!
        function delete(obj)
            if(obj.Handle ~= -1)
                calllib(Splinter.getInstance().get_alias(), obj.Destructor_function, obj.Handle);
            end
        end
    end
end

