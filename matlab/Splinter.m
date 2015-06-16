classdef Splinter < handle
    properties (Access = protected)
        Alias
    end

    methods(Access = private)
        function obj = Splinter
            obj.Alias = 'libSplinter';
        end 
    end
    
    methods(Static)
        function singleObj = getInstance
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = Splinter;
            end
            singleObj = localObj;
        end
    end
    
    methods
        function load(obj, libFile, headerFile, alias)
            if(obj.is_loaded())
                fprintf('Splinter is already loaded as %s\n', obj.Alias);
                fprintf('If you wish to change the alias you should unload and then load the library with the new alias.\n');
            else
                if(exist('alias', 'var'))
                   obj.Alias = alias;
                end
                
                fprintf('Loading Splinter as %s\n', obj.Alias);
                loadlibrary(libFile, headerFile, 'alias', obj.Alias)
            end
        end
        
        function r = call(obj, func, varargin)
            if(~obj.is_loaded())
                error('You need to load the library with Splinter.getInstance().load(libFile, headerFile [, alias]) before making calls to it!\n'); 
            end
            
            % strcmp == 1 if identical
            if(strcmp(func, ''))
                error('This function has not been implemented yet!'); 
            end
            
            r = calllib(obj.Alias, func, varargin{:});
            
            % Check the internal error flag (set when trying to reference
            % an object that does not exist, which may happen when creating
            % an object, reloading the library, then trying to use the
            % object. The internal reference will then be gone, and the
            % reference invalid).
            if(calllib(obj.Alias, 'get_error'))
               error('Internal library error: Did you use an object belonging to a previous load of the library?'); 
            end
        end
        
        function unload(obj)
           unloadlibrary(obj.Alias);
        end
        
        function r = is_loaded(obj)
           r = libisloaded(obj.Alias);
        end
        
        function alias = get_alias(obj)
           alias = obj.Alias; 
        end
    end
end

