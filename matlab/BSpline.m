classdef BSpline
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = private)
        Handle
    end
    
    methods
        % Constructor. Creates an instance of the BSpline class in the
        % library, using the samples in dataTable.
        % type is an integer specifying the degree of the bspline
        % type can be any one of: 1, 3, 4
        function obj = BSpline(dataTable, type)
           obj.Handle = calllib('splinter', 'bspline_init', dataTable.get_handle(), type);
        end
        
        % Evaluate the b-spline at x
        % x should be a real number (not array) within the domain of the b-spline
        function r = eval(obj, x)
           r = calllib('splinter', 'bspline_eval', obj.Handle, x);
        end
        
        % Destructor. Deletes the internal BSpline object.
        % Note that calling any of the functions after this,
        % WILL make MatLab crash!
        function delete(obj)
           calllib('splinter', 'bspline_delete', obj.Handle); 
        end
    end
end

