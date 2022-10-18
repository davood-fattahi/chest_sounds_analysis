classdef cvpartition2
    properties(GetAccess = public, SetAccess = public, Dependent=false)
       %TYPE  Validation partition type. 
       %   The Type property is a string with value  'kfold', 'holdout',
       %   'leaveout' or 'resubstitution'. It indicates the type of
       %   validation partition.
        Type; 
        
       %NUMTESTSETS Number of test sets.
       %   The NumTestSets property is an integer indicating the number of
       %   folds in K-fold and Leave-one-out. The value is one in Holdout
       %   and Resubstitution.
        NumTestSets;
        
        %TRAINSIZE Size of each training set.
        %   The TrainSize property indicates the size of each training set.
        %   It is a numeric vector in K-fold and Leave-one-out, or a
        %   numeric scalar in Holdout and Resubstitution.
        TrainSize;
        
        %TESTSIZE Size of each test set.
        %   The TestSize property indicates the size of each test set. It
        %   is a numeric vector in K-fold and Leave-one-out, or a numeric
        %   scalar in Holdout and Resubstitution.
        TestSize;
        
        %NUMOBSERVATIONS Number of observations.
        %   The NumObservations property is a numeric scalar holding the number of
        %   observations, including observations with missing GROUP values
        %   if GROUP is provided.
         NumObservations;       
         
        full_test;
        full_train;
        
        Impl; 
    end

    methods
        function n = get.Type(this)
            n = this.Type;
        end
        function n = get.NumTestSets(this)
            n = this.NumTestSets;
        end
        function n = get.TrainSize(this)
            n = this.TrainSize;
        end
        function n = get.TestSize(this)
            n = this.TestSize;
        end
        function n = get.NumObservations(this)
            n = this.NumObservations;
        end        
       
        function trainIndices = training(this,I)
        trainIndices = this.full_train{I}; 
        end

        function testIndices = test(this,I)
        testIndices = this.full_test{I}; 
        end
        
        function disp(this)
            %disp(this)
        end
        
        function obj = cvpartition2(Type,NumTestSets,TrainSize,TestSize,NumObservations,full_test,full_train)
            if nargin > 0
                obj.Type = Type;
                obj.NumTestSets = NumTestSets;
                obj.TrainSize = TrainSize; 
                obj.TestSize = TestSize; 
                obj.NumObservations = NumObservations; 
                obj.full_test = full_test;
                obj.full_train = full_train; 
            end
        end
    end % public methods block

    
    methods(Hidden = true)
        
        % Methods that we inherit, but do not want
        function a = fields(varargin),     throwUndefinedError(); end
        function a = ctranspose(varargin), throwUndefinedError(); end
        function a = transpose(varargin),  throwUndefinedError(); end
        function a = permute(varargin),    throwUndefinedError(); end
        function a = reshape(varargin),    throwUndefinedError(); end
        function a = cat(varargin),        throwNoCatError(); end
        function a = horzcat(varargin),    throwNoCatError(); end
        function a = vertcat(varargin),    throwNoCatError(); end
    end
    methods(Hidden = true, Static = true)
        function a = empty(varargin)
            error(message('stats:cvpartition:NoEmptyAllowed', upper( mfilename )));
        end
    end
    
    methods(Hidden, Static, Access='public')
         function cv = loadobj(cv)
             currentVersion = 1; % Should be in sync with the cv.Version property
             if isstruct(cv)
                 % The incompatibility introduced in 16b will always reach
                 % this point with a struct (not a object). cvpartition
                 % constructor has a backdoor for this case, so we do not
                 % need to take any other preventive action now.
                 cv.BckCmpBackdoorConstructor = true;
                 cv = cvpartition(cv);
             elseif isempty(cv.Version) || (cv.Version < currentVersion)
                 % This branch should never be reached
                 error(message('stats:cvpartition:FwrdCompatibility'))
             end
         end
    end % hidden static public
   
end % classdef

function throwNoCatError()
error(message('stats:cvpartition:NoCatAllowed', upper( mfilename )));
end

function throwUndefinedError()
st = dbstack;
name = regexp(st(2).name,'\.','split');
error(message('stats:cvpartition:UndefinedFunction', name{ 2 }, mfilename));
end
