classdef DataTable < handle
    properties (Access = private)
        Handle
        
        Samples
        X_dim
        Num_samples
        
        Is_synced
    end
    
    methods
        % Constructor. Creates an instance of the DataTable class in the
        % library.
        function obj = DataTable()
            obj.Handle = Splinter.getInstance().call('datatable_init');
            obj.Samples = [];
            obj.Num_samples = 0;
            obj.Is_synced = true;
        end
        
        % Add a sample at x with value y.
        % x should be a real number (not array)!
        % Multi-dimensional support is in the library but not
        % the MatLab interface (yet).
        function r = add_sample(obj, x, y)
            r = Splinter.getInstance().call('datatable_add_sample', obj.Handle, x, length(x), y);
        end
        
        % Reserve memory for the Samples. It is *HIGHLY* recommended that
        % you do this if you know how many samples you are going to insert,
        % as there are significant time savings to be had by not
        % reallocating the matrix every time a new sample is inserted.
        % You may also benefit from doing this even if you guess the number
        % of samples, it does not have to be exact.
        function preallocate(obj, numRows, numCols)
            obj.Samples = NaN(numRows, numCols);
        end
        
        % Add a sample to the Samples matrix in this class. They are sent
        % to the library by the finish-method. This is to avoid a call via
        % calllib every time a sample is inserted. This method is about
        % twice as fast as inserting a sample one by one into the library.
        function add_samples(obj, x, y)
            obj.Samples(obj.Num_samples+1,:) = [x y];
            if(obj.Num_samples == 0)
                obj.X_dim = length(x);
            end
            obj.Num_samples = obj.Num_samples + 1;
            obj.Is_synced = false;
        end
        
        % Send the samples to the library
        function finish(obj)
            % Make a libpointer to avoid copying the matrix
            if(obj.Is_synced == false)
                samplePtr = libpointer('doublePtr', obj.Samples);
                Splinter.getInstance().call('datatable_add_samples', obj.Handle, samplePtr, obj.Num_samples, obj.X_dim) 
                obj.Is_synced = true;
            end
        end
        
        % Internal use only
        function r = get_handle(obj)
            % Sync with library in case the DataTable is about to get used
            %obj.finish()
            r = obj.Handle;
        end
        
        % Destructor. Deletes the internal DataTable object.
        function delete(obj)
            Splinter.getInstance().call('datatable_delete', obj.Handle) 
        end
    end
end

