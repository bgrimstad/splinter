% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef PSpline < Approximant
    properties (Access = protected)
        Handle
        
        Constructor_function = 'pspline_init';
        Constructor_load_function = 'pspline_load_init';
    end

    methods
        % Constructor. Creates an instance of the PSpline class in the
        % library, using the samples in dataTable.
        % lambda is the smoothing parameter, usually a small number
        % default: 0.03
        function obj = PSpline(dataTableOrFilename, lambda)
            if(~exist('lambda', 'var'))
                lambda = 0.03;
            end
            
            % Set to -1 so we don't try to delete the library instance in case type is invalid
            obj.Handle = -1;
            
            if(ischar(dataTableOrFilename))
                filename = dataTableOrFilename;
                
                obj.Handle = Splinter.getInstance().call(obj.Constructor_load_function, filename);
            else
                dataTable = dataTableOrFilename;
                
                obj.Handle = Splinter.getInstance().call(obj.Constructor_function, dataTable.get_handle(), lambda);
            end
        end
    end
end