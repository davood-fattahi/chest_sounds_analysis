function make(opt)
fprintf('=> Building SVOREX.\n');
if nargin < 1
    try
      % This part is for OCTAVE
      if (exist ('OCTAVE_VERSION', 'builtin'))
		% Use -std=c++11 for newer versions of Octave
        if ispc
          setenv('CFLAGS','-std=gnu99 -O3')
          setenv('CC','gcc')
        else
          setenv('CFLAGS','-O3 -fstack-protector-strong -Wformat -Werror=format-security')
        end
        mex -O3 mainSvorex.c alphas.c cachelist.c datalist.c def_settings.c kcv.c loadfile.c ordinal_takestep.c setandfi.c smo_kernel.c smo_routine.c smo_settings.c smo_timer.c svc_predict.c -output svorex
        delete *.o
      % This part is for MATLAB
      else
          if ispc
            mex COMPFLAGS="\$COMPFLAGS -std=c99 -O3" -largeArrayDims mainSvorex.c alphas.c cachelist.c datalist.c def_settings.c kcv.c loadfile.c ordinal_takestep.c setandfi.c smo_kernel.c smo_routine.c smo_settings.c smo_timer.c svc_predict.c -output svorex
          else
            mex CFLAGS="\$CFLAGS -std=c99 -O3" -largeArrayDims mainSvorex.c alphas.c cachelist.c datalist.c def_settings.c kcv.c loadfile.c ordinal_takestep.c setandfi.c smo_kernel.c smo_routine.c smo_settings.c smo_timer.c svc_predict.c -output svorex
          end
      end
    catch err
        fprintf('Error: %s failed (line %d)\n', err.stack(1).file, err.stack(1).line);
        disp(err.message);
        fprintf('=> Please check README for detailed instructions.\n');
    end
elseif nargin == 1
    switch lower(opt)
        case 'clean'
            delete *.o
        case 'cleanall'
            delete *.o
            delete *.mexa64
        otherwise
            error('make option "%s" not recognized', opt)
    end
end
