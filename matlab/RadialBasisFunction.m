% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef RadialBasisFunction < Approximant
    properties (Access = protected)
        Handle
        
        Constructor_function = 'rbf_init';
    end

    methods
        function obj = RadialBasisFunction(dataTable, type, normalized)
            % Set to -1 so we don't try to delete the library instance in case type is invalid
            obj.Handle = -1;
            
            % Default to Thin plate spline
            if(~exist('type', 'var'))
                type_index = 1;
            else
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
            
            obj.Handle = Splinter.getInstance().call(obj.Constructor_function, dataTable.get_handle(), type_index,  normalized);
        end
    end
end

