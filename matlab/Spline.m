% loadlibrary('C:/Users/Anders/Documents/GitHub/splinter/Debug/splinter-matlab-1-2.dll', 'C:/Users/Anders/Documents/GitHub/splinter/splinter/include/matlab.h', 'alias', 'splinter')
classdef (Abstract = true) Spline < handle
    properties (Constant)
        Splinter_alias = 'splinter'
    end
    
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
            r = calllib(obj.Splinter_alias, obj.Eval_function, obj.Handle, x, length(x));
        end
        
        % Evaluate the Jacobian at x
        % x should be a real number (not array) within the domain of the spline
        function r = evalJacobian(obj, x)
            temp = calllib(obj.Splinter_alias, obj.EvalJacobian_function, obj.Handle, x, length(x));
            reshape(temp, 1, length(x))
            r = temp.value;
        end
        
        % Evaluate the Hessian at x
        % x should be a real number (not array) within the domain of the spline
        function r = evalHessian(obj, x)
            temp = calllib(obj.Splinter_alias, obj.EvalHessian_function, obj.Handle, x, length(x));
            reshape(temp, length(x), length(x))
            r = temp.value;
        end
        
        function save(obj, fileName)
            calllib(obj.Splinter_alias, obj.Save_function, obj.Handle, fileName);
        end

        % Destructor. Deletes the internal BSpline object.
        % Note that calling any of the methods after this,
        % WILL make MatLab crash!
        function delete(obj)
            if(obj.Handle ~= -1)
                calllib(Spline.Splinter_alias, obj.Destructor_function, obj.Handle);
            end
        end
    end
end

