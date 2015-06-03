classdef BSpline < Spline
    properties (Access = protected)
        Handle
        
        Constructor_function = 'bspline_init';
        Destructor_function = 'bspline_delete';
        Eval_function = 'bspline_eval';
        EvalJacobian_function = 'bspline_eval_jacobian';
        EvalHessian_function = 'bspline_eval_hessian';
        Save_function = 'bspline_save';
        Load_function = 'bspline_load';
    end

    methods
        % Constructor. Creates an instance of the BSpline class in the
        % library, using the samples in dataTable.
        % type is an enumeration constant (from BSplineType)
        % specifying the degree of the b-spline.
        function obj = BSpline(dataTable, type)
            % Set to -1 so we don't try to delete the library instance in case type is invalid
            obj.Handle = -1;
            
            % These values are somewhat arbitrary, the important thing is
            % that the MatLab front end agrees with the C back end
            type_index = -1;
            switch type
                case BSplineType.Linear
                    type_index = 0;
                case BSplineType.Quadratic
                    type_index = 1;
                case BSplineType.Quadratic_free
                    type_index = 2;
                case BSplineType.Cubic
                    type_index = 3;
                case BSplineType.Cubic_free
                    type_index = 4;
            end
            
            if(type_index == -1)
                error('type should be an enumeration constant of type BSplineType!')
            else
                obj.Handle = calllib(obj.Splinter_alias, obj.Constructor_function, dataTable.get_handle(), type_index);
                
                if(obj.Handle == 0)
                    error('Could not create BSpline!');
                end
            end
        end
    end
end

