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
            fprintf('Calling %s with %d args\n', func, length(nargin));
            if(~obj.is_loaded())
                error('You need to load the library with Splinter.getInstance().load(libFile, headerFile [, alias]) before making calls to it!\n'); 
            end
            
            r = calllib(obj.Alias, func, varargin{:});
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

