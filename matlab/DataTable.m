classdef DataTable < handle
    properties (Constant)
       Splinter_alias = 'splinter'; 
    end
    
    properties (Access = private)
        Handle
    end
    
    methods
        % Constructor. Creates an instance of the DataTable class in the
        % library.
        function obj = DataTable()
           obj.Handle = calllib(obj.Splinter_alias, 'datatable_init');
        end
        
        % Add a sample at x with value y.
        % x should be a real number (not array)!
        % Multi-dimensional support is in the library but not
        % the MatLab interface (yet).
        function r = add_sample(obj, x, y)
           r = calllib(obj.Splinter_alias, 'datatable_add_sample', obj.Handle, x, length(x), y);
        end
        
        % Destructor. Deletes the internal DataTable object.
        function delete(obj)
           calllib(obj.Splinter_alias, 'datatable_delete', obj.Handle) 
        end
        
        % Internal use only
        function r = get_handle(obj)
            r = obj.Handle;
        end
    end
end

