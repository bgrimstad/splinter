% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef PolynomialRegression < Approximant
    properties (Access = protected)
        Handle
        
        Constructor_function = 'polynomial_regression_init';
        Constructor_load_function = 'polynomial_regression_load_init';
    end

    methods
        % Constructor. Creates an instance of the BSpline class in the
        % library, using the samples in dataTable.
        % type is an enumeration constant (from BSplineType)
        % specifying the degree of the b-spline.
        function obj = PolynomialRegression(dataTableOrFilename, degree)
            % Set to -1 so we don't try to delete the library instance in case type is invalid
            obj.Handle = -1;
            
            if(ischar(dataTableOrFilename))
                filename = dataTableOrFilename;
                
                obj.Handle = Splinter.getInstance().call(obj.Constructor_load_function, filename);
            else
                dataTable = dataTableOrFilename;
                
                degree_vec = degree;
                if(isscalar(degree))
                    degree_vec = ones(1,dataTable.get_num_variables()) * degree;
                end
                
                obj.Handle = Splinter.getInstance().call(obj.Constructor_function, dataTable.get_handle(), degree_vec, length(degree_vec));
            end
        end
    end
end

