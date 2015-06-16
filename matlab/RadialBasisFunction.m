classdef RadialBasisFunction < Approximant
    properties (Access = protected)
        Handle
        
        Constructor_function = 'rbf_init';
        Destructor_function = 'rbf_delete';
        Eval_function = 'rbf_eval';
        EvalJacobian_function = '';
        EvalHessian_function = '';
        Save_function = 'rbf_save';
        Load_function = 'rbf_load';
    end

    methods
        function obj = RadialBasisFunction(dataTable, type, normalized)
            % Set to -1 so we don't try to delete the library instance in case type is invalid
            obj.Handle = -1;
            
            % Default to Thin plate spline
            type_index = 1;
            switch type
                case RadialBasisFunctionType.Thin_plate_spline
                    type_index = 1;
                case RadialBasisFunctionType.Multiquadric
                    type_index = 2;
                case RadialBasisFunctionType.Inverse_quadric
                    type_index = 3;
                case RadialBasisFunctionType.Inverse_multiquadric
                    type_index = 4;
                case RadialBasisFunctionType.Gaussian
                    type_index = 4;
            end
            
            if(~exist('normalized', 'var'))
               normalized = 0;
            end
            
            % Make sure that normalized is 0 or 1
            if(normalized ~= 0)
                normalized = 1;
            else
                normalized = 0;
            end
            
            obj.Handle = calllib(Splinter.getInstance().get_alias(), obj.Constructor_function, dataTable.get_handle(), type_index,  normalized);

            if(obj.Handle == 0)
                error('Could not create RadialBasisFunction!');
            end
        end
    end
end

