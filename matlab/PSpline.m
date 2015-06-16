classdef PSpline < Approximant
    properties (Access = protected)
        Handle
        
        Constructor_function = 'pspline_init';
        Destructor_function = 'pspline_delete';
        Eval_function = 'pspline_eval';
        EvalJacobian_function = 'pspline_eval_jacobian';
        EvalHessian_function = 'pspline_eval_hessian';
        Save_function = 'pspline_save';
        Load_function = 'pspline_load';
    end

    methods
        % Constructor. Creates an instance of the PSpline class in the
        % library, using the samples in dataTable.
        % lambda is the smoothing parameter, usually a small number
        % default: 0.03
        function obj = PSpline(dataTable, lambda)
            if(~exist('lambda', 'var'))
                lambda = 0.03;
            end
            
            % Set to -1 so we don't try to delete the library instance in case type is invalid
            obj.Handle = -1;
            
            obj.Handle = Splinter.getInstance().call(obj.Constructor_function, dataTable.get_handle(), lambda);
        end
    end
end

