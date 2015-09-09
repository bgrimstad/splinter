% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef DataTable < handle
    properties (Access = private)
        Handle
        
        X_dim
        
        % Samples store the samples that are awaiting transfer to the back
        % end, Num_samples is tracking how many samples we are storing
        Samples
        Num_samples
    end
    
    methods
        % Constructor. Creates an instance of the DataTable class in the
        % library.
        function obj = DataTable(filename)
            obj.Samples = [];
            obj.Num_samples = 0;
            obj.X_dim = 0;
            
            if(exist('filename', 'var') && ischar(filename))
                obj.Handle = Splinter.getInstance().call('datatable_load_init', filename);
                obj.X_dim = obj.get_num_variables();
            else
                obj.Handle = Splinter.getInstance().call('datatable_init');
            end
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
        
        % Add a sample to the DataTable.
        function add_sample(obj, x, y)
            if(obj.X_dim == 0)
                obj.X_dim = length(x);
            end
            
            if(length(x) ~= obj.X_dim)
                error('Dimension of new sample is inconsistent with previous samples!')
            end
            
            obj.Samples(obj.Num_samples+1,:) = [x y];
            obj.Num_samples = obj.Num_samples + 1;
        end
        
        function r = get_num_variables(obj)
            obj.finish()
            r = Splinter.getInstance().call('datatable_get_num_variables', obj.Handle); 
        end
        
        function r = get_num_samples(obj)
            obj.finish()
            r = Splinter.getInstance().call('datatable_get_num_samples', obj.Handle); 
        end
        
        % Send the samples to the library
        function finish(obj)
            % Make a libpointer to avoid copying the matrix
            if(obj.Num_samples > 0)
                samplePtr = libpointer('doublePtr', obj.Samples);
                
                % We need to pass the number of rows in the matrix to the
                % library because of the way MatLab stores matrices in
                % memory. It stores them column-major, which means that we
                % need to know the number of rows to reconstruct the data
                % in the back end. The number of rows does not necessarily
                % equal the number of samples, because we can add some
                % data, send it to the library, then add some more. If the
                % number of samples in the last sequence is less than that
                % in the first, the number of rows will be more than the
                % number of samples, and the library will get trash data if
                % it uses the number of samples as if it equaled the number
                % of rows.
                temp = size(obj.Samples);
                numRows = temp(1);
                Splinter.getInstance().call('datatable_add_samples', obj.Handle, samplePtr, obj.Num_samples, obj.X_dim, numRows)
                
                % This number is the number of samples we are storing
                % temporarily in the MatLab front end, and that are
                % awaiting transfer to the back end. Therefore we need to
                % set it to 0 here.
                obj.Num_samples = 0;
            end
        end
        
        function save(obj, filename)
           Splinter.getInstance().call('datatable_save', obj.Handle, filename)
        end
        
        % Internal use only
        function r = get_handle(obj)
            % Sync with library in case the DataTable is about to get used
            obj.finish()
            r = obj.Handle;
        end
        
        % Destructor. Deletes the internal DataTable object.
        function delete(obj)
            if(obj.Handle ~= -1)
                Splinter.getInstance().call('datatable_delete', obj.Handle) 
            end
        end
    end
end