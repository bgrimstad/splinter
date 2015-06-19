% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef BSpline < Approximant
    properties (Access = protected)
        Handle
        
        Constructor_function = 'bspline_init';
        Destructor_function = 'bspline_delete';
        Eval_function = 'bspline_eval';
        EvalJacobian_function = 'bspline_eval_jacobian';
        EvalHessian_function = 'bspline_eval_hessian';
        GetNumVariables_function = 'bspline_get_num_variables';
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
                    type_index = 1;
                case BSplineType.Quadratic
                    type_index = 2;
                case BSplineType.Cubic
                    type_index = 3;
                case BSplineType.Quartic
                    type_index = 4;
            end
            
            if(type_index == -1)
                error('type should be an enumeration constant of type BSplineType!')
            else
                obj.Handle = Splinter.getInstance().call(obj.Constructor_function, dataTable.get_handle(), type_index);
            end
        end
    end
end

