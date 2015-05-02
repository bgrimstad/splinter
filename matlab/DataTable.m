classdef DataTable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        Handle
    end
    
    methods
        % Constructor. Creates an instance of the DataTable class in the
        % library.
        function obj = DataTable()
           obj.Handle = calllib('splinter', 'datatable_init');
        end
        
        % Add a sample at x with value y.
        % x should be a real number (not array)!
        % Multi-dimensional support is in the library but not
        % the MatLab interface (yet).
        function r = add_sample(obj, x, y)
           r = calllib('splinter', 'datatable_add_sample', obj.Handle, x, y);
        end
        
        % Destructor. Deletes the internal DataTable object.
        function delete(obj)
           calllib('splinter', 'datatable_delete', obj.Handle) 
        end
        
        % Internal use only
        function r = get_handle(obj)
            r = obj.Handle;
        end
    end
end

